#!/usr/bin/env bash
# Fetch the DNA of each PRDM9 locus found by an earlier run of run_toga_prdm9.sh,
# WITHOUT downloading the TOGA2 annotation files again.
#
# Uses what is already on disk:
#   prdm9_toga_raw/<short>.PRDM9.loci.tsv   (the loci found last time)
#   genome_2bit_urls.tsv                    (built by make_genome_urls.py)
# and writes the same files run_toga_prdm9.sh would have written:
#   prdm9_locus_seq/<short>__<scaffold>_<start>_<end>.fa
#
#   bash fetch_locus_seq.sh            # all assemblies that have a loci file
#   bash fetch_locus_seq.sh myoMyo6    # just these ones
set -uo pipefail

GENOME_URLS="genome_2bit_urls.tsv"
RAWDIR="./prdm9_toga_raw"
SEQDIR="./prdm9_locus_seq"
TMPDIR="./prdm9_seq_tmp"
LOG="prdm9_seq.log"
SEQ_FLANK=100000                 # keep equal to SEQ_FLANK in run_toga_prdm9.sh

mkdir -p "$SEQDIR" "$TMPDIR"; touch "$LOG"
cleanup() { rm -rf "$TMPDIR"; }
trap cleanup EXIT INT TERM

command -v twoBitToFa >/dev/null || { echo "ERROR: twoBitToFa not found" >&2; exit 1; }
[[ -s "$GENOME_URLS" ]] || { echo "ERROR: $GENOME_URLS missing - run: python3 make_genome_urls.py" >&2; exit 1; }

if [[ $# -gt 0 ]]; then
    files=(); for s in "$@"; do files+=("$RAWDIR/${s}.PRDM9.loci.tsv"); done
else
    shopt -s nullglob; files=("$RAWDIR"/*.PRDM9.loci.tsv)
fi
[[ ${#files[@]} -gt 0 ]] || { echo "No *.PRDM9.loci.tsv in $RAWDIR" >&2; exit 1; }

n_ok=0; n_skip=0; n_fail=0
for lf in "${files[@]}"; do
    [[ -s "$lf" ]] || { echo "no loci file: $lf" >&2; continue; }
    species=$(basename "$lf" .PRDM9.loci.tsv)
    url=$(awk -F'\t' -v s="$species" '$1==s{print $2; exit}' "$GENOME_URLS")
    if [[ -z "$url" ]]; then
        echo -e "${species}\tNO_GENOME_URL\t-\tnot listed in ${GENOME_URLS}" >> "$LOG"
        echo "  ${species}: no genome URL" >&2; n_fail=$((n_fail+1)); continue
    fi

    # scaffold lengths: chrom.sizes sits next to the 2bit file
    sizes="${TMPDIR}/${species}.sizes"
    curl -sSf --retry 2 "$(dirname "$url")/chrom.sizes" -o "$sizes" 2>/dev/null || : > "$sizes"

    udc="${TMPDIR}/udc_${species}"; mkdir -p "$udc"
    while IFS=$'\t' read -r chr s e _; do
        [[ -z "${chr:-}" ]] && continue
        size=$(awk -v c="$chr" '$1==c{print $2; exit}' "$sizes")
        ss=$(( s > SEQ_FLANK ? s - SEQ_FLANK : 0 )); se=$(( e + SEQ_FLANK ))
        [[ -n "$size" && "$se" -gt "$size" ]] && se=$size
        fa="${SEQDIR}/${species}__${chr}_${ss}_${se}.fa"
        if [[ -s "$fa" ]]; then n_skip=$((n_skip+1)); continue; fi
        err="${TMPDIR}/err"
        if twoBitToFa -udcDir="$udc" -seq="$chr" -start="$ss" -end="$se" "$url" "$fa" 2>"$err" \
           && [[ -s "$fa" ]]; then
            sed -i "1s/.*/>${species}|${chr}|${ss}/" "$fa"
            echo -e "${species}\tSEQ_OK\t${chr}:${ss}-${se}" >> "$LOG"
            n_ok=$((n_ok+1))
        else
            rm -f "$fa"
            echo -e "${species}\tSEQ_FAILED\t${chr}:${ss}-${se}\t$(tr '\n' ' ' < "$err" | cut -c1-140)" >> "$LOG"
            n_fail=$((n_fail+1))
        fi
    done < "$lf"
    rm -rf "$udc"
    echo "  ${species}: done" >&2
done

echo "fetched ${n_ok}, already had ${n_skip}, failed ${n_fail}. Details in ${LOG}" >&2
[[ $n_fail -gt 0 ]] && echo "reasons:" >&2 && grep SEQ_FAILED "$LOG" | cut -f4 | sort | uniq -c >&2
exit 0
