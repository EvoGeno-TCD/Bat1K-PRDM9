#!/usr/bin/env bash
# Pull PRDM9 projections from TOGA2 bigBeds, one assembly at a time.
#
# PRDM9's zinc-finger array resembles hundreds of other zinc-finger genes, so
# TOGA reports many PRDM9 projections across the genome. ALL of them are kept
# in prdm9_toga_raw/<short>.PRDM9.bed (for the status/loss summary), but only
# the ORTHOLOGOUS loci (see ORTHOLOGUE FILTER below) get:
#   - their neighbouring genes (+/-NEIGH_FLANK)  -> synteny check
#   - their DNA sequence (+/-SEQ_FLANK)          -> re-annotation with miniprot
# The downloaded bigBed is cached only while one assembly is processed, then deleted.
set -uo pipefail

URLS="bat_toga_urls.tsv"                 # short<TAB>url
PARSER="parse_toga_prdm9.py"             # included; standard library only
GENE_REGEX='prdm9'                       # case-insensitive gene name match

# ---- ORTHOLOGUE FILTER ----------------------------------------------------
MIN_CHAIN=0.5     # keep projections with chain orthology score (column 16) >= this
MAX_LOCI=5        # at most this many loci per assembly (highest chain score first)
# Projections whose name carries an extra flag (e.g. ...#PRDM9#123#retro) are never
# used as loci. If nothing passes, the single best-scoring projection is used and
# the log says FILTER_FALLBACK.

NEIGH_FLANK=1000000                      # 1 Mb either side (clipped to scaffold ends)
SEQ_FLANK=100000                         # 100 kb either side
FETCH_SEQ=1                              # 0 = skip sequence download
# Sequence source 1 (preferred): the Hiller lab's own genome 2bit files, on the
# same server as the annotations. Build the list first with:
#     python3 make_genome_urls.py
GENOME_URLS="genome_2bit_urls.tsv"       # short<TAB>https://.../<assembly>.2bit
# Sequence source 2 (fallback, used only when 1 is unavailable for an assembly):
# NCBI, which works only when the scaffold names are NCBI accessions (CM..., NC_...).
USE_NCBI_FALLBACK=1
NCBI_API_KEY="${NCBI_API_KEY:-}"         # optional; raises NCBI's rate limit
EFETCH_PAUSE=0.4                         # seconds between NCBI requests (3/s limit without key)

OUTDIR="./prdm9_toga_clean"
RAWDIR="./prdm9_toga_raw"
SEQDIR="./prdm9_locus_seq"
TMPDIR="./prdm9_tmp"
LOG="prdm9_extract.log"

mkdir -p "$OUTDIR" "$RAWDIR" "$SEQDIR" "$TMPDIR"
touch "$LOG"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT INT TERM

[[ -f "$PARSER" ]] || { echo "ERROR: $PARSER not found in $(pwd)" >&2; exit 1; }
if [[ "$FETCH_SEQ" == 1 && ! -s "$GENOME_URLS" ]]; then
    echo "WARNING: $GENOME_URLS not found - locus sequences will be tried from NCBI only," >&2
    echo "         which fails whenever scaffold names are not NCBI accessions." >&2
    echo "         Build it first:  python3 make_genome_urls.py" >&2
fi

purge_cache() {   # $1 = assembly dir (HLxxx), $2 = species short name
    rm -rf "${TMPDIR:?}/udc_${2}"
    rm -rf "/tmp/udcCache/http/genome.senckenberg.de/data/${1}" \
           "/tmp/udcCache/https/genome.senckenberg.de/data/${1}"
}

fetch_ncbi() {   # $1 seqname $2 start(0-based) $3 end -> stdout FASTA
    local key=""; [[ -n "$NCBI_API_KEY" ]] && key="&api_key=${NCBI_API_KEY}"
    curl -sSf --retry 3 --retry-delay 2 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${1}&rettype=fasta&strand=1&seq_start=$(( $2 + 1 ))&seq_stop=${3}${key}"
    sleep "$EFETCH_PAUSE"
}

while IFS=$'\t ' read -r species url _; do
    [[ -z "${species:-}" || "$species" =~ ^# || "$species" == short ]] && continue
    species="${species//$'\r'/}"; url="${url//$'\r'/}"
    asm_dir=$(basename "$(dirname "$url")")
    out_tsv="${OUTDIR}/${species}.PRDM9.clean.tsv"
    raw_bed="${RAWDIR}/${species}.PRDM9.bed"
    nbr_bed="${RAWDIR}/${species}.PRDM9.neighbours.bed"
    loci_out="${RAWDIR}/${species}.PRDM9.loci.tsv"
    [[ -s "$out_tsv" || -e "${RAWDIR}/${species}.none" ]] && { echo "Skipping ${species}" >&2; continue; }

    echo "Processing ${species}" >&2
    tmp_bed="${TMPDIR}/${species}.bed"; tmp_tsv="${TMPDIR}/${species}.parsed.tsv"
    udc="${TMPDIR}/udc_${species}"; mkdir -p "$udc"
    sizes="${TMPDIR}/${species}.sizes"

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
    bigBedInfo -udcDir="$udc" -chroms "$url" 2>/dev/null \
        | awk 'NF==3 && $3 ~ /^[0-9]+$/ {print $1"\t"$3}' > "$sizes"

    # ---- orthologous loci: no extra flag in the name, chain score >= MIN_CHAIN ----
    awk -F'\t' -v m="$MIN_CHAIN" '{n=split($4,p,"#")} n<=3 && ($16+0)>=m' "$tmp_bed" \
        > "${TMPDIR}/${species}.orth.bed"
    if [[ ! -s "${TMPDIR}/${species}.orth.bed" ]]; then
        sort -t$'\t' -k16,16gr "$tmp_bed" | head -1 > "${TMPDIR}/${species}.orth.bed"
        echo -e "${species}\tFILTER_FALLBACK\tno projection passed MIN_CHAIN=${MIN_CHAIN}; used best one" >> "$LOG"
    fi
    # merge overlapping projections into loci, keep best chain score per locus, top MAX_LOCI
    sort -t$'\t' -k1,1 -k2,2n "${TMPDIR}/${species}.orth.bed" | awk 'BEGIN{FS=OFS="\t"}
        $1!=c || $2>e { if (c!="") print c,s,e,b; c=$1; s=$2; e=$3; b=$16+0; next }
        { if ($3>e) e=$3; if ($16+0>b) b=$16+0 }
        END{ print c,s,e,b }' \
      | sort -t$'\t' -k4,4gr | head -n "$MAX_LOCI" > "$loci_out"
    n_all=$(wc -l < "$tmp_bed"); n_loci=$(wc -l < "$loci_out")
    echo "  ${n_all} PRDM9 projections; ${n_loci} orthologous locus/loci kept" >&2

    # ---- neighbours: compact columns only, windows clipped to the scaffold ----
    : > "$nbr_bed"
    while IFS=$'\t' read -r chr s e _; do
        size=$(awk -v c="$chr" '$1==c{print $2}' "$sizes")
        ws=$(( s > NEIGH_FLANK ? s - NEIGH_FLANK : 0 )); we=$(( e + NEIGH_FLANK ))
        [[ -n "$size" && "$we" -gt "$size" ]] && we=$size
        bigBedToBed -udcDir="$udc" -chrom="$chr" -start="$ws" -end="$we" "$url" stdout \
            | awk -v L="${chr}:${s}-${e}" 'BEGIN{FS=OFS="\t"}
                  {print L,$1,$2,$3,$4,$6,$16,$23}' >> "$nbr_bed" \
            || echo -e "${species}\tNEIGHBOURS_FAILED\t${chr}:${s}-${e}" >> "$LOG"
    done < "$loci_out"
    purge_cache "$asm_dir" "$species"

    # ---- locus sequence: Senckenberg 2bit first, NCBI fallback -------------------
    if [[ "$FETCH_SEQ" == 1 ]]; then
        tb_url=""
        [[ -s "$GENOME_URLS" ]] && tb_url=$(awk -F'\t' -v s="$species" '$1==s{print $2; exit}' "$GENOME_URLS")
        while IFS=$'\t' read -r chr s e _; do
            size=$(awk -v c="$chr" '$1==c{print $2}' "$sizes")
            ss=$(( s > SEQ_FLANK ? s - SEQ_FLANK : 0 )); se=$(( e + SEQ_FLANK ))
            [[ -n "$size" && "$se" -gt "$size" ]] && se=$size
            fa="${SEQDIR}/${species}__${chr}_${ss}_${se}.fa"
            [[ -s "$fa" ]] && continue
            ok=0; why=""
            err="${TMPDIR}/${species}.seqerr"
            if [[ -n "$tb_url" ]]; then
                # -seq/-start/-end, not 'file:seq:start-end': the URL contains colons
                if twoBitToFa -udcDir="${TMPDIR}/udc_seq" -seq="$chr" -start="$ss" -end="$se" \
                        "$tb_url" "$fa" 2>"$err" && [[ -s "$fa" ]]; then
                    ok=1
                else
                    why="2bit: $(tr '\n' ' ' < "$err" | cut -c1-120)"
                fi
                rm -rf "${TMPDIR}/udc_seq"
            else
                why="no entry in ${GENOME_URLS}"
            fi
            if [[ $ok == 0 && "$USE_NCBI_FALLBACK" == 1 ]]; then
                if fetch_ncbi "$chr" "$ss" "$se" > "$fa" 2>"$err" \
                   && head -c1 "$fa" | grep -q '>' \
                   && [[ $(grep -v '>' "$fa" | tr -d '\n' | wc -c) -gt 1000 ]]; then
                    ok=1
                else
                    why="${why}; ncbi: $(tr '\n' ' ' < "$err" | cut -c1-120)"
                fi
            fi
            if [[ $ok == 1 ]]; then
                sed -i "1s/.*/>${species}|${chr}|${ss}/" "$fa"
            else
                rm -f "$fa"
                echo -e "${species}\tSEQ_FAILED\t${chr}:${ss}-${se}\t${why}" >> "$LOG"
            fi
            rm -f "$err"
        done < "$loci_out"
    fi

    # ---- parse all PRDM9 rows ----------------------------------------------------
    if ! python3 "$PARSER" "$tmp_bed" "$tmp_tsv"; then
        echo -e "${species}\tPARSER_FAILED\t${url}" >> "$LOG"; continue
    fi
    awk -v sp="$species" 'BEGIN{FS=OFS="\t"} NR==1{$1="species"; print; next} {$1=sp; print}' \
        "$tmp_tsv" > "$out_tsv"
    rm -f "$tmp_bed" "$tmp_tsv" "${TMPDIR}/${species}".*
    echo -e "${species}\tOK\tprojections=${n_all}\tloci=${n_loci}" >> "$LOG"
done < "$URLS"

echo "Done. Summary in ${LOG}:" >&2
cut -f2 "$LOG" | sort | uniq -c >&2
