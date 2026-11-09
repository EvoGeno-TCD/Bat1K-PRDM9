# PRDM9 in bats — pipeline

## Run order
```bash
python3 make_assembly_metadata.py toga_table.txt bat_assemblies.tsv   # tech / haplotype / individual
bash run_toga_prdm9.sh        # TOGA PRDM9 rows + flanking genes + locus sequence; UDC cache purged per assembly
bash reannotate_prdm9.sh      # miniprot + read-through ZF exon  (needs query_prdm9.faa)
bash prdm9_downstream.sh      # everything below
```

## Before the first run, check
| What | Where |
|---|---|
| Status column number in the bigBed (`bigBedInfo -as <url>`) | `STATUS_COL` in `prdm9_downstream.sh` |
| 2bit URL pattern for locus sequence (test one with `twoBitInfo`) | `TWOBIT_PATTERN` in `run_toga_prdm9.sh` |
| The 3 non-HL assemblies (myoBra1, myoDav1, pteAle1) resolve | `bat_toga_urls.tsv` |
| `tech` = `unknown` for 57 assemblies — fill in where you know | `bat_assemblies.tsv` |
| Query proteins: human PRDM9 (Q9NQV7) + 1–3 bat PRDM9s | `query_prdm9.faa` |
| Published dated tree + synonym map | `BIG_TREE`, `NAME_MAP` |
| Predicted motifs + genome paths (erosion test) | `motifs.tsv`, `genomes.tsv` |

## Scripts
| Script | Analysis |
|---|---|
| `make_assembly_metadata.py` | sequencing tech (heuristic), haplotype, individual (ToLID) per assembly |
| `run_toga_prdm9.sh` | extraction; neighbours ±1 Mb; locus ±100 kb; cache deletion |
| `reannotate_prdm9.py/.sh` | ZF exon read from genome to first stop; repairs fake introns; gap flags |
| `domain_table.py` | KRAB / SSXRD / SET / C2H2 by accession; domain order |
| `prdm9_zf_analysis.py` | ZF array, counts, ZF0, contact residues (−1,+2,+3,+6), unit FASTAs, non-ZF CDS |
| `prdm9_qc.py` | status incl. losses; synteny score; PRDM9 gate; hap1 vs hap2; version agreement; one representative/species |
| `zf_dnds.py repeats` | pooled Nei–Gojobori dN/dS among repeats, contact vs non-contact, bootstrap CIs |
| `zf_dnds.py pairs` | pairwise dN/dS on codon alignment; dS saturation flag |
| `concerted_test.py` | nearest-neighbour + within/between permutation tests |
| `prune_tree.py` | prune published tree to representatives, relabel for HyPhy |
| `zf_count_asr.R` | ZF-count ancestral states, phylogenetic signal, branch changes > allelic noise (untested) |
| `motif_erosion.py` | own-genome motif depletion vs relatives (Markov-1 expectation) |
| `fasta_select.py` | subset/rename FASTAs to representatives |

## Key outputs (prdm9_analysis/)
`qc.status_by_species.tsv` (losses) · `qc.master.tsv` · `qc.haplotypes.tsv` · `qc.representatives.tsv` ·
`prdm9.per_zf.tsv` · `zf.repeats.tsv` · `nonzf.pairs.tsv` · `concerted.summary.txt` ·
`nonzf.{busted,meme,absrel}.json` · `zfcount.branch_changes.tsv` · `erosion.erosion.tsv`
