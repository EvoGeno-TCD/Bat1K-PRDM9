#!/usr/bin/env bash
# Re-annotate PRDM9 at every locus fetched by run_toga_prdm9.sh.
# Needs: miniprot (>=0.12), python3.
#
# QUERY: protein FASTA of PRDM9 queries. Use human PRDM9 (UniProt Q9NQV7)
# plus 1-3 good bat PRDM9 proteins from long-read assemblies (e.g. from a
# first pass) - a bat query aligns much better than human across the KRAB/SET.
set -euo pipefail

QUERY="query_prdm9.faa"
SEQDIR="./prdm9_locus_seq"
GFFDIR="./prdm9_miniprot_gff"
OUT="./prdm9_reannot"
THREADS=4
mkdir -p "$GFFDIR"

[[ -s "$QUERY" ]] || { echo "Put PRDM9 query proteins in $QUERY" >&2; exit 1; }

for fa in "$SEQDIR"/*.fa; do
    base=$(basename "$fa" .fa)
    gff="$GFFDIR/${base}.gff"
    [[ -s "$gff" ]] && continue
    # -N 3: keep a few alternative models; -p 0.5: report weaker hits too
    miniprot -t "$THREADS" --gff -N 3 -p 0.5 "$fa" "$QUERY" > "$gff" 2>/dev/null \
        || echo "miniprot failed: $base" >&2
done

python3 reannotate_prdm9.py --gff-dir "$GFFDIR" --seq-dir "$SEQDIR" -o "$OUT"
echo "Flags:" >&2; cut -f18 "${OUT}.report.tsv" | tail -n +2 | tr ',' '\n' | sort | uniq -c >&2
# Next: prdm9_downstream.sh uses ${OUT}.faa / ${OUT}.fna as its candidates.
