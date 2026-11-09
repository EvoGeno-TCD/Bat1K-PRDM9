#!/usr/bin/env bash
# Pull PRDM9 projections from TOGA2 bigBeds, one assembly at a time.
# For each assembly, while the remote file is still cached, also grab:
#   - every TOGA projection within +/-NEIGH_FLANK of each PRDM9 hit (synteny check)
#   - the genomic sequence of the PRDM9 locus +/-SEQ_FLANK (for re-annotation)
# Then delete the UDC cache and keep only the small outputs.
set -uo pipefail

URLS="bat_toga_urls.tsv"                 # short<TAB>url
PARSER="parse_toga_ifne.py"              # your existing parser
GENE_REGEX='prdm9'                       # case-insensitive
NEIGH_FLANK=1000000                      # 1 Mb either side
SEQ_FLANK=100000                         # 100 kb either side
FETCH_SEQ=1                              # 0 to skip sequence download
# 2bit location for each assembly - VERIFY ONE BY HAND FIRST, e.g.
#   twoBitInfo http://genome.senckenberg.de/data/HLmyoMyo6/HLmyoMyo6.2bit stdout | head
TWOBIT_PATTERN='http://genome.senckenberg.de/data/{ASM}/{ASM}.2bit'

OUTDIR="./prdm9_toga_clean"
RAWDIR="./prdm9_toga_raw"
SEQDIR="./prdm9_locus_seq"
TMPDIR="./prdm9_tmp"
UDC_ROOT="${TMPDIR}/udc"
LOG="prdm9_extract.log"

mkdir -p "$OUTDIR" "$RAWDIR" "$SEQDIR" "$TMPDIR"
touch "$LOG"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT

purge_cache() {   # $1 = assembly dir (HLxxx), $2 = species short name
    rm -rf "${UDC_ROOT:?}/${2}"
    rm -rf "/tmp/udcCache/http/genome.senckenberg.de/data/${1}"
}

while IFS=$'\t ' read -r species url _; do
    [[ -z "${species:-}" || "$species" =~ ^# || "$species" == short ]] && continue
    species="${species//$'\r'/}"; url="${url//$'\r'/}"
    asm_dir=$(basename "$(dirname "$url")")
    out_tsv="${OUTDIR}/${species}.PRDM9.clean.tsv"
    raw_bed="${RAWDIR}/${species}.PRDM9.bed"
    nbr_bed="${RAWDIR}/${species}.PRDM9.neighbours.bed"
    [[ -s "$out_tsv" || -e "${RAWDIR}/${species}.none" ]] && { echo "Skipping ${species}" >&2; continue; }

    echo "Processing ${species}" >&2
    tmp_bed="${TMPDIR}/${species}.bed"; tmp_tsv="${TMPDIR}/${species}.parsed.tsv"
    udc="${UDC_ROOT}/${species}"; mkdir -p "$udc"

    if ! bigBedToBed -udcDir="$udc" "$url" stdout \
         | awk -v re="#${GENE_REGEX}([#.]|\$)" 'BEGIN{FS=OFS="\t"} tolower($4) ~ re' \
         > "$tmp_bed"; then
        echo -e "${species}\tDOWNLOAD_FAILED\t${url}" >> "$LOG"
        purge_cache "$asm_dir" "$species"; continue
    fi

    if [[ ! -s "$tmp_bed" ]]; then
        echo -e "${species}\tNO_PRDM9\t${url}" >> "$LOG"
        touch "${RAWDIR}/${species}.none"
        purge_cache "$asm_dir" "$species"; continue
    fi
    cp "$tmp_bed" "$raw_bed"

    # ---- neighbourhood (cheap: file already cached) ------------------------
    : > "$nbr_bed"
    sort -k1,1 -k2,2n "$tmp_bed" | awk 'BEGIN{FS=OFS="\t"}
        $1!=c || $2>e {if (c!="") print c,s,e; c=$1; s=$2; e=$3; next}
        {if ($3>e) e=$3} END{print c,s,e}' > "${TMPDIR}/${species}.loci"
    while IFS=$'\t' read -r chr s e; do
        ws=$(( s > NEIGH_FLANK ? s - NEIGH_FLANK : 0 )); we=$(( e + NEIGH_FLANK ))
        bigBedToBed -udcDir="$udc" -chrom="$chr" -start="$ws" -end="$we" "$url" stdout \
            | awk -v L="${chr}:${s}-${e}" 'BEGIN{FS=OFS="\t"}{print L,$0}' >> "$nbr_bed" \
            || echo -e "${species}\tNEIGHBOURS_FAILED\t${chr}:${s}-${e}" >> "$LOG"
    done < "${TMPDIR}/${species}.loci"
    purge_cache "$asm_dir" "$species"

    # ---- locus sequence from the assembly's 2bit ---------------------------
    if [[ "$FETCH_SEQ" == 1 ]]; then
        tb="${TWOBIT_PATTERN//\{ASM\}/$asm_dir}"
        while IFS=$'\t' read -r chr s e; do
            ss=$(( s > SEQ_FLANK ? s - SEQ_FLANK : 0 )); se=$(( e + SEQ_FLANK ))
            fa="${SEQDIR}/${species}__${chr}_${ss}_${se}.fa"
            if twoBitToFa -udcDir="$udc" "${tb}:${chr}:${ss}-${se}" "$fa" 2>/dev/null \
               && [[ -s "$fa" ]]; then
                # twoBitToFa clips at contig end; record the real offset in the header
                sed -i "1s/.*/>${species}|${chr}|${ss}/" "$fa"
            else
                rm -f "$fa"
                echo -e "${species}\tSEQ_FAILED\t${chr}:${ss}-${se}" >> "$LOG"
            fi
        done < "${TMPDIR}/${species}.loci"
        rm -rf "${UDC_ROOT:?}/${species}"
        rm -rf "/tmp/udcCache/http/genome.senckenberg.de/data/${asm_dir}"
    fi

    # ---- parse -------------------------------------------------------------
    if ! python "$PARSER" "$tmp_bed" "$tmp_tsv"; then
        echo -e "${species}\tPARSER_FAILED\t${url}" >> "$LOG"; continue
    fi
    awk -v sp="$species" 'BEGIN{FS=OFS="\t"} NR==1{$1="species"; print; next} {$1=sp; print}' \
        "$tmp_tsv" > "$out_tsv"
    rm -f "$tmp_bed" "$tmp_tsv"
    echo -e "${species}\tOK\t$(($(wc -l < "$out_tsv") - 1))" >> "$LOG"
done < "$URLS"

echo "Done. Summary in ${LOG}:" >&2
cut -f2 "$LOG" | sort | uniq -c >&2
