#!/usr/bin/env bash
# PRDM9 bat analysis driver. Order of the whole pipeline:
#   0. python3 make_assembly_metadata.py <toga_table.txt> bat_assemblies.tsv
#   1. bash run_toga_prdm9.sh          (TOGA rows, neighbours, locus sequence)
#   2. bash reannotate_prdm9.sh        (miniprot + read-through ZF exon)
#   3. bash prdm9_downstream.sh        (this script)
# Tools: interproscan.sh, mafft, pal2nal.pl, iqtree2, hyphy, Rscript (ape, phytools)
# Steps whose inputs are missing are skipped with a message, not failed.
set -euo pipefail

META="bat_assemblies.tsv"
RAWDIR="./prdm9_toga_raw"
CLEANDIR="./prdm9_toga_clean"        # parser output (inactivating_mutations)
MIN_CHAIN=0.5                        # keep equal to MIN_CHAIN in run_toga_prdm9.sh
ZF_EXON="auto"                       # reference exon number of the ZF exon, or auto
REANNOT="./prdm9_reannot"            # .faa .fna .report.tsv from reannotate_prdm9.sh
STATUS_COL=23                        # check: bigBedInfo -as <url>
BIG_TREE=""                          # published dated tree (Newick/NEXUS), e.g. Upham 2019 MCC
NAME_MAP=""                          # optional species<TAB>tree_label synonyms
MOTIFS="motifs.tsv"                  # species<TAB>IUPAC motif (from a ZF predictor)
GENOMES="genomes.tsv"                # species<TAB>genome.fa (for erosion test)
W="./prdm9_analysis"
THREADS=8
mkdir -p "$W"
have() { command -v "$1" >/dev/null 2>&1; }
skip() { echo "SKIP: $*" >&2; }

FAA="${REANNOT}.faa"; FNA="${REANNOT}.fna"
[[ -s "$FAA" ]] || { echo "Run reannotate_prdm9.sh first ($FAA missing)" >&2; exit 1; }

# 1. Domain architecture ------------------------------------------------------
if [[ ! -s "$W/prdm9.interpro.tsv" ]] && have interproscan.sh; then
    # InterProScan refuses any sequence containing '*' (a stop codon). Copies with
    # an internal stop are pseudogenes we still want domains for, so replace '*'
    # with 'X' (the IUPAC "unknown residue"). Every position stays where it was,
    # so the coordinates in prdm9.per_zf.tsv still line up.
    IPR_IN="$W/prdm9_for_interpro.faa"
    awk '/^>/ {print; next}
         {n += gsub(/\*/, "X"); print}
         END {if (n) printf("replaced %d stop codon(s) with X for InterProScan\n", n) > "/dev/stderr"}' \
        "$FAA" > "$IPR_IN"
    interproscan.sh -i "$IPR_IN" -f tsv -cpu "$THREADS" \
        -appl Pfam,SMART,ProSiteProfiles,ProSitePatterns -o "$W/prdm9.interpro.tsv"
fi
DOM_ARG=()
if [[ -s "$W/prdm9.interpro.tsv" ]]; then
    python3 domain_table.py "$W/prdm9.interpro.tsv" "$FAA" > "$W/prdm9.domains.tsv"
    DOM_ARG=(--domains "$W/prdm9.domains.tsv")
else skip "InterProScan not available - PRDM9 gate will not check KRAB/SSXRD/SET"; fi

# 2. ZF arrays ---------------------------------------------------------------
python3 prdm9_zf_analysis.py "$FAA" --cds "$FNA" -o "$W/prdm9"

# 3. QC: status/losses, synteny, gate, haplotypes, representatives -------------
python3 prdm9_qc.py --meta "$META" --rawdir "$RAWDIR" --status-col "$STATUS_COL" \
    --cleandir "$CLEANDIR" --zf-exon "$ZF_EXON" --min-chain-score "$MIN_CHAIN" \
    --zf "$W/prdm9.per_protein.tsv" "${DOM_ARG[@]}" \
    --reannot "${REANNOT}.report.tsv" -o "$W/qc"
REPS="$W/qc.representatives.tsv"
# representatives that pass the gate (used for all cross-species analyses)
awk -F'\t' 'NR==1 || $7=="PRDM9"' "$REPS" > "$W/reps.pass.tsv"
echo "species with a gated PRDM9 representative: $(($(wc -l < "$W/reps.pass.tsv") - 1))" >&2

# 4. Non-ZF domains: alignment, codon alignment, saturation, gene tree ----------
python3 fasta_select.py "$W/prdm9.nonzf.faa" "$W/reps.pass.tsv" --rename > "$W/nonzf.reps.faa"
python3 fasta_select.py "$W/prdm9.nonzf.fna" "$W/reps.pass.tsv" --rename > "$W/nonzf.reps.fna"
if have mafft; then
    mafft --localpair --maxiterate 1000 --thread "$THREADS" "$W/nonzf.reps.faa" > "$W/nonzf.aln.faa"
    if have pal2nal.pl; then
        pal2nal.pl "$W/nonzf.aln.faa" "$W/nonzf.reps.fna" -output fasta > "$W/nonzf.codon.fa"
        python3 zf_dnds.py pairs "$W/nonzf.codon.fa" -o "$W/nonzf" || skip "pairwise dN/dS"
        echo "saturated pairs (dS>1.5): $(grep -c SATURATED "$W/nonzf.pairs.tsv" || true)" >&2
    else skip "pal2nal.pl"; fi
    have iqtree2 && iqtree2 -s "$W/nonzf.aln.faa" -m MFP -B 1000 -T AUTO --prefix "$W/nonzf" -redo
else skip "mafft"; fi

# 5. ZF repeats: dN/dS contact vs non-contact (pooled within arrays) -----------
python3 fasta_select.py "$W/prdm9.zf_units.fna" "$W/reps.pass.tsv" --units > "$W/zf_units.reps.fna"
python3 zf_dnds.py repeats "$W/zf_units.reps.fna" -o "$W/zf" --boot 200 \
    || skip "repeat dN/dS failed"

# 6. Concerted evolution ----------------------------------------------------------
python3 concerted_test.py "$W/prdm9.zf_units.faa" --keep "$W/reps.pass.tsv" \
    --meta "$META" --perm 1000 -o "$W/concerted" || skip "concerted test (needs >=2 arrays)"
if have mafft && have iqtree2; then
    python3 fasta_select.py "$W/prdm9.zf_units.fna" "$W/reps.pass.tsv" --units --rename \
        > "$W/zf_units.reps.named.fna"
    mafft --globalpair --maxiterate 1000 "$W/zf_units.reps.named.fna" > "$W/zf_units.aln.fna"
    iqtree2 -s "$W/zf_units.aln.fna" -m MFP -B 1000 -T AUTO --prefix "$W/zf_units" -redo
fi

# 7. Species tree, tree comparison, selection, ZF-count mapping ----------------
if [[ -n "$BIG_TREE" && -s "$BIG_TREE" ]]; then
    python3 prune_tree.py "$BIG_TREE" "$W/reps.pass.tsv" ${NAME_MAP:+--map "$NAME_MAP"} -o "$W/bats"
    [[ -s "$W/bats.missing.txt" ]] && echo "unmatched species: see $W/bats.missing.txt (add to NAME_MAP)" >&2
    TREE="$W/bats.short.nwk"

    # gene tree vs species tree (Robinson-Foulds)
    if have iqtree2 && [[ -s "$W/nonzf.treefile" ]]; then
        iqtree2 -rf "$TREE" "$W/nonzf.treefile" --prefix "$W/rf_nonzf_vs_species" -redo || true
    fi

    # selection on non-ZF domains: recombination first, then gene/site/branch tests
    if have hyphy && [[ -s "$W/nonzf.codon.fa" ]]; then
        hyphy gard   --alignment "$W/nonzf.codon.fa" --output "$W/nonzf.gard.json" || true
        # If GARD finds breakpoints, rerun the tests below on $W/nonzf.codon.fa.best-gard
        hyphy busted --alignment "$W/nonzf.codon.fa" --tree "$TREE" --output "$W/nonzf.busted.json"
        hyphy meme   --alignment "$W/nonzf.codon.fa" --tree "$TREE" --output "$W/nonzf.meme.json"
        hyphy absrel --alignment "$W/nonzf.codon.fa" --tree "$TREE" --output "$W/nonzf.absrel.json"
    else skip "HyPhy or codon alignment"; fi

    if have Rscript; then
        Rscript zf_count_asr.R "$TREE" "$W/reps.pass.tsv" "$W/qc.master.tsv" \
            "$W/qc.haplotypes.tsv" "$W/zfcount" || skip "R ASR failed (check packages)"
    fi
else skip "species tree (set BIG_TREE)"; fi

# 8. Hotspot erosion ------------------------------------------------------------
# Predict motifs from $W/prdm9.array.faa with a C2H2 recognition-code predictor,
# write $MOTIFS, list genomes in $GENOMES, then:
if [[ -s "$MOTIFS" && -s "$GENOMES" ]]; then
    python3 motif_erosion.py "$MOTIFS" "$GENOMES" --meta "$META" --same-family -o "$W/erosion" \
        || skip "erosion test failed"
else skip "erosion test (needs $MOTIFS and $GENOMES)"; fi

echo "Done. Start with: $W/qc.status_by_species.tsv  $W/qc.master.tsv  $W/qc.haplotypes.tsv" >&2
