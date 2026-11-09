# PRDM9 evolution in bats

Pipeline that finds the gene **PRDM9** in 160 bat
genome assemblies (108 species), determines how reliable each copy is, and then
tests how PRDM9 has evolved across bats.

Read sections 1–3 before running anything. They explain what the pipeline does
and why. Section 4 sets up the software, sections 5-6 cover running it, and
section 7 explains how to read the results.

---

## Contents
1. [Background: why PRDM9?](#1-background-why-prdm9)
2. [The questions we are asking](#2-the-questions-we-are-asking)
3. [How the pipeline works (overview)](#3-how-the-pipeline-works-overview)
4. [Software you need](#4-software-you-need)
5. [Files you must prepare](#5-files-you-must-prepare)
6. [Running it, step by step](#6-running-it-step-by-step)
7. [What the outputs mean](#7-what-the-outputs-mean)
8. [Common problems](#8-common-problems)
9. [Known limitations](#9-known-limitations)
10. [Glossary](#10-glossary)
11. [Key references](#11-key-references)

---

## 1. Background: why PRDM9?

During meiosis, chromosomes swap pieces (**recombination**). This doesn't
happen at random places. In humans, mice and many other mammals, it is
concentrated in short stretches called **hotspots**, and **PRDM9** decides
where those hotspots are.

PRDM9 is a protein with several parts (**domains**), arranged in this order
from start to end:

| Domain | What it does |
|---|---|
| **KRAB** | Interacts with other proteins |
| **SSXRD** | Short repression domain; function less clear |
| **SET** (also called PR/SET) | Adds chemical marks (H3K4me3, H3K36me3) to histones, which flags the site for recombination |
| **C2H2 zinc-finger (ZF) array** | A row of repeated "fingers". Each finger grips ~3 DNA bases, so together they decide which DNA **motif** PRDM9 binds |

Three facts make PRDM9 interesting to study in evolution:

- **The ZF array evolves very fast.** It is encoded by a tandem repeat
  (~84 bp per finger) in a single exon, the last one. Fingers are gained,
  lost and copied between repeats.
- **Only a few amino acids in each finger touch DNA.** These are at helix
  positions **−1, +2, +3 and +6**. Changes here change the binding motif.
  Earlier work found these positions under strong positive selection
  (Oliver et al. 2009).
- **Hotspots destroy themselves.** Recombination at a hotspot tends to
  replace the bound motif with the partner chromosome's version, so over time
  the motif is lost from the genome. PRDM9 then has to evolve to bind new
  motifs. This is the "Red Queen" idea, and it predicts a species' own motif
  should be **under-represented** in its own genome.

PRDM9 has also been **lost or partly lost** several times in vertebrates
(dogs, birds and others; Baker et al. 2017). Having all of KRAB + SSXRD + SET
intact predicts whether PRDM9 still controls hotspots. So whether each bat has
a complete PRDM9 is a result in its own right.

---

## 2. The questions we are asking

| # | Question | Main output |
|---|---|---|
| 1 | Does each bat have a complete PRDM9 (KRAB, SSXRD, SET, ZF array)? Any losses? | `qc.status_by_species.tsv`, `prdm9.domains.tsv` |
| 2 | How many zinc fingers does each species have, and where on the bat tree did arrays expand or shrink? | `prdm9.per_protein.tsv`, `zfcount.branch_changes.tsv` |
| 3 | How different are the DNA-contact residues between species (i.e. do they bind different motifs)? | `prdm9.per_zf.tsv`, `prdm9.binding_codes.tsv` |
| 4 | Is there positive selection, and is it stronger at DNA-contact residues? | `zf.repeats.tsv`, HyPhy `.json` files |
| 5 | Are fingers within one species more similar to each other than to other species' fingers (**concerted evolution**)? | `concerted.summary.txt`, `concerted.pairs.tsv` |
| 6 | Is each species' predicted motif depleted in its own genome (**hotspot erosion**)? | `erosion.erosion.tsv` |

---

## 3. How the pipeline works (overview)

```
 TOGA2 annotations (online)            your inputs
 ─────────────────────────             ───────────
        │                              query_prdm9.faa, published bat tree,
        ▼                              motifs.tsv, genomes.tsv
 [Step 0] make_assembly_metadata.py  → which assemblies are long-read, which are haplotypes
        │
 [Step 1] run_toga_prdm9.sh          → PRDM9 rows from TOGA, nearby genes, DNA of the locus
        │                              (big download per genome, deleted straight after)
 [Step 2] reannotate_prdm9.sh        → rebuild the PRDM9 protein from the DNA itself
        │
 [Step 3] prdm9_downstream.sh        → domains, ZF arrays, QC, selection, concerted
                                       evolution, tree mapping, erosion
```

### Why each step exists

**Step 0: assembly metadata.** Genome assemblies differ in quality. Older
assemblies built from **short reads** often collapse repeats, and the ZF array
is a repeat. So a short-read assembly can report too few fingers. We record each
assembly's technology so that ZF counts are only trusted from **long-read**
assemblies. Some species have two assemblies of the same animal, one per
chromosome copy (**hap1/hap2**). These are two alleles, not replicates, and
PRDM9 is very variable even within species.

**Step 1: extraction from TOGA2.** TOGA2 (Hiller lab) has already annotated
thousands of genes in these genomes by projecting human/mouse genes onto them.
We pull out:
- every PRDM9 projection and its TOGA **status** (Intact, Partially intact, Lost …);
- every gene within 1 Mb either side, so we can check that the copy sits in the
  expected genomic neighbourhood (**synteny**). This separates the real gene
  from copies elsewhere in the genome;
- the DNA sequence of the PRDM9 region ±100 kb.

Each annotation file is downloaded into a temporary cache, which is
**deleted after every assembly** to save disk space.

**Step 2: re-annotation from DNA.** TOGA works well for most of the gene but is
weakest exactly where we care most: the repetitive ZF exon. It can drop or merge
fingers. So we align a known PRDM9 protein to the DNA of each region (with
**miniprot**). We then read the final exon directly from the DNA until the stop
codon. That way a bat with more fingers than the human query still gets all of
them. Assembly gaps (runs of N) in that exon are flagged.

**Step 3: analysis.** See section 7 for each output.
1. **Domains:** InterProScan finds KRAB, SSXRD, SET and C2H2 domains.
2. **ZF arrays:** finds every finger, separates the C-terminal tandem array from
   the lone upstream finger (**ZF0**), and pulls out the contact residues.
3. **QC:** combines TOGA status, mutations, synteny, domains and assembly
   quality. A copy is only called "PRDM9" if it has SET + SSXRD, has a ZF array,
   sits in the expected neighbourhood and has a clean ZF exon. It then picks
   **one representative assembly per species** for cross-species analyses
   (avoids counting one species several times).
4. **Non-ZF domains:** aligned across species (MAFFT), converted to a codon
   alignment (pal2nal) and used for selection tests and a gene tree.
5. **dN/dS among repeats:** compares fingers within each array, separately for
   contact and non-contact codons.
6. **Concerted evolution test.**
7. **Species tree work:** prune a published bat tree to our species, compare it
   with the PRDM9 gene tree, run HyPhy selection tests, and map ZF counts onto it.
8. **Hotspot erosion test.**

### Why ZF repeats are not aligned across species for selection tests
Standard selection tests (codeml, HyPhy) assume each alignment column is the
**same position** in every species. Finger number 5 in one bat is not
necessarily related to finger 5 in another, because fingers are gained, lost and
copied. So we:
- run the standard cross-species tests only on the **non-ZF** part of the gene;
- test the ZF array by comparing fingers **within** each array (step 5),
  following Oliver et al. (2009).

---

## 4. Softwares you will need

The easiest route is a conda/mamba environment with the bioconda channel:

```bash
mamba create -n prdm9 -c conda-forge -c bioconda \
    python=3.11 numpy beautifulsoup4 \
    ucsc-bigbedtobed ucsc-bigbedinfo ucsc-twobittofa ucsc-twobitinfo \
    miniprot mafft pal2nal iqtree hyphy \
    r-base r-ape r-phytools \
    ncbi-datasets-cli
mamba activate prdm9
```

| Tool | Used in | Needed for | Required? |
|---|---|---|---|
| bash, awk, coreutils | all `.sh` scripts | glue | **yes** |
| Python ≥ 3.8 | all `.py` scripts | everything | **yes** |
| numpy | `concerted_test.py`, `motif_erosion.py` | statistics | **yes** |
| beautifulsoup4 | `parse_toga_ifne.py` (lab parser) | parsing TOGA rows | **yes** |
| `bigBedToBed` (UCSC) | step 1 | reading TOGA annotation files over the web | **yes** |
| `twoBitToFa` (UCSC) | step 1 | downloading the PRDM9 region | **yes** |
| `bigBedInfo`, `twoBitInfo` (UCSC) | checks only | checking settings before the first run | recommended |
| miniprot ≥ 0.12 | step 2 | protein-to-genome alignment | **yes** |
| InterProScan 5 | step 3.1 | finding domains | strongly recommended* |
| MAFFT | step 3.4, 3.6 | alignments | recommended |
| pal2nal.pl | step 3.4 | codon alignment | recommended |
| IQ-TREE 2 (`iqtree2`) | steps 3.4, 3.6, 3.7 | gene trees, tree comparison | recommended |
| HyPhy | step 3.7 | selection tests | recommended |
| R + `ape` + `phytools` | step 3.7 | mapping ZF counts on the tree | recommended |
| NCBI `datasets` | before step 3.8 | downloading whole genomes | only for erosion test |

\* **InterProScan** is large (tens of GB with its databases) and needs Java. The
bioconda version sometimes lacks databases, so the official install from the
InterPro website is safer. Without it the pipeline still runs, but the PRDM9
check can't look at domains.

**IQ-TREE version note:** the scripts call `iqtree2`. If you only have IQ-TREE 3
(`iqtree3`), either install version 2 or make an alias:
`ln -s "$(which iqtree3)" ~/bin/iqtree2`.

**Steps that skip themselves:** if a recommended tool is missing, step 3 prints
`SKIP: <reason>` and carries on. Read the end of the log to see what didn't run.

**Disk and network:** Step 1 downloads one large annotation file per assembly,
which can be several GB. It needs internet access to `genome.senckenberg.de`.
The cache is deleted after each assembly, so you only need room for one at a
time. The erosion test will need whole genomes (~2 GB each).

---

## 5. Files to prepare

Put **all scripts and these files in one folder**, and run everything from
that folder.

### Already provided

| File | What it is |
|---|---|
| `bat_toga_urls.tsv` | One line per bat assembly: short name, TAB, TOGA2 annotation URL. Used by step 1. |
| `bat_assemblies.tsv` | Made by step 0: species, family, technology, haplotype, individual for each assembly. |

### You must add

| File | Format | Used by | How to get it |
|---|---|---|---|
| Published dated bat/mammal **species** tree | Newick or NEXUS | step 3.7 (`BIG_TREE=` in `prdm9_downstream.sh`) | One summary tree (e.g. the MCC tree) is enough. Tip labels must start with `Genus_species`. This is a tree of the **species**, not of the gene; the pipeline builds the gene trees itself. |
| `name_map.tsv` (optional) | `Genus species<TAB>exact tree tip label` | step 3.7 (`NAME_MAP=`) | Only for species whose names changed (e.g. *Neoeptesicus*, *Cnephaeus*, *Dermanura*). After the first run, look at `prdm9_analysis/bats.missing.txt` to see which ones need it. |
| `motifs.tsv` | `short_name<TAB>IUPAC motif` | step 3.8 (`MOTIFS=`) | Run a C2H2 zinc-finger binding predictor (e.g. the Persikov & Singh predictor) on `prdm9_analysis/prdm9.array.faa` and write one motif per species. Use the **short** assembly names (e.g. `myoNatt2A`). |
| `genomes.tsv` | `short_name<TAB>/path/to/genome.fa(.gz)` | step 3.8 (`GENOMES=`) | See below. |

**Making `genomes.tsv`:**
- **Labels:** must be **identical** to those in `motifs.tsv`.
- **Which genomes:** list each representative genome (the one the motif came
  from), plus at least one other genome from the same bat family. Motifs are
  only compared within a family.
- **Downloading:** for assemblies with an accession, use
  `datasets download genome accession GCA_xxxxxxx.x --include genome`.
  For assemblies without an accession, convert the TOGA2 hub's genome file:
  `twoBitToFa http://genome.senckenberg.de/data/HLxxx/HLxxx.2bit xxx.fa`.
  Then delete `/tmp/udcCache/...` afterwards.
- **Starting the file from the representatives:**
  ```bash
  awk -F'\t' 'NR==FNR{acc[$1]=$6; next} FNR>1{print $4"\t/data/genomes/"acc[$4]".fna.gz"}' \
      bat_assemblies.tsv prdm9_analysis/reps.pass.tsv > genomes.tsv
  ```

### Settings to check before the first run

| Setting | In file | What it is | How to check |
|---|---|---|---|
| `STATUS_COL=23` | `prdm9_downstream.sh` | Which column of the TOGA rows holds the status | `bigBedInfo -as <any URL from bat_toga_urls.tsv>`. The lab parser also reads status from this column. |
| `TWOBIT_PATTERN` | `run_toga_prdm9.sh` | Web address of each genome's sequence file | `twoBitInfo http://genome.senckenberg.de/data/HLmyoMyo6/HLmyoMyo6.2bit stdout \| head` should list scaffolds |
| `ZF_EXON="auto"` | `prdm9_downstream.sh` | Which reference exon holds the ZF array | Step 3 prints the auto-detected number; check it against the human PRDM9 transcript |
| `STATUS_RANK` list | top of `prdm9_qc.py` | Expected TOGA status words | After step 1: `cut -f15 prdm9_toga_clean/*.tsv \| sort \| uniq -c` |
| `tech` column | `bat_assemblies.tsv` | Sequencing technology | 57 assemblies are `unknown`. Where you can find out (NCBI assembly page), set one of `long_read_phased`, `long_read`, `short_read_HiC`, `short_read` (exact spelling). |
| `THREADS` | `reannotate_prdm9.sh`, `prdm9_downstream.sh` | CPUs to use | Match your machine |

---

## 6. Running the pipeline - step by step

Run each step in order and check its output before moving on. On a shared
server, run long steps inside `screen` or `tmux`, or as a cluster job.

### Step 0: assembly metadata (seconds)
```bash
python3 make_assembly_metadata.py toga_species_table.txt bat_assemblies.tsv
```
`toga_species_table.txt` is the TOGA2 species list (the tab-separated table of
species, assembly, taxonomy ID and lineage). `bat_assemblies.tsv` is already
provided, so you only need to re-run this if that table changes.

### Step 1: extract from TOGA2 (hours; network-bound)
```bash
bash run_toga_prdm9.sh
cut -f2 prdm9_extract.log | sort | uniq -c       # summary
```
Codes in `prdm9_extract.log`:

| Code | Meaning | What to do |
|---|---|---|
| `OK` | PRDM9 found and parsed | nothing |
| `NO_PRDM9` | TOGA has no PRDM9 for this assembly | a possible loss; noted in QC |
| `DOWNLOAD_FAILED` | URL didn't work | check the URL (likely for myoBra1, myoDav1, pteAle1) |
| `SEQ_FAILED` | region sequence not downloaded | check `TWOBIT_PATTERN` |
| `NEIGHBOURS_FAILED` | nearby genes not fetched | synteny check missing for that copy |
| `PARSER_FAILED` | lab parser crashed | usually beautifulsoup4 missing |

The script **resumes**: it skips any assembly that already has
`prdm9_toga_clean/<name>.PRDM9.clean.tsv` or `prdm9_toga_raw/<name>.none`. To
redo one, delete those files first.

### Step 2: re-annotate from DNA (minutes)
```bash
bash reannotate_prdm9.sh
```
At the end it prints a count of flags (see section 7, `prdm9_reannot.report.tsv`).
When you add bat proteins to `query_prdm9.faa`, delete `prdm9_miniprot_gff/`
before re-running, or old alignments are reused.

### Step 3: analysis (minutes to hours, depending on tools)
First set `BIG_TREE` (and `NAME_MAP`, if needed) at the top of
`prdm9_downstream.sh`, then:
```bash
bash prdm9_downstream.sh 2>&1 | tee downstream.log
grep SKIP downstream.log           # which parts did not run, and why
```
Everything is written to `prdm9_analysis/`.
`prdm9_analysis/prdm9.interpro.tsv` is reused if present; delete it after
changing step 2.

### Re-running only part of step 3
Every script can be run on its own; `python3 <script>.py -h` explains its
options. For example:
```bash
python3 concerted_test.py prdm9_analysis/prdm9.zf_units.faa \
    --keep prdm9_analysis/reps.pass.tsv --meta bat_assemblies.tsv -o test_run
```

---

## 7. What the outputs mean

All paths are inside `prdm9_analysis/` unless stated otherwise.
**Suggested reading order:** 1 → 2 → 3 → 4 → then the analysis outputs.

### 1. `qc.status_by_species.tsv`: does each species have PRDM9?
One row per species.

| Column | Meaning |
|---|---|
| `best_status_longread` / `best_status_other` | Best TOGA status among long-read / other assemblies |
| `evidence_assembly` | Assembly the evidence comes from (long-read preferred) |
| `loss_evidence` | Type of damage in that assembly's best PRDM9 copy (table below) |
| `interpretation` | Plain-language reading of `loss_evidence` |
| `note` | e.g. "candidate loss/truncation" |

`loss_evidence`, from strongest to weakest evidence of a real loss:

| Value | Meaning |
|---|---|
| `outside_zf_exon` | Frameshift/stop/splice mutation in a non-ZF exon. **Strongest evidence of real loss.** |
| `deleted_exons` | Exon(s) deleted outside the ZF exon. Loss evidence, but check. |
| `missing_exons_only` | Exons not found, often just assembly gaps. Weak. |
| `zf_exon_only` | Problems only in the ZF exon. Often a TOGA artefact; compare with the step 2 protein. |
| `non_inactivating_only` / `none` | No sign of loss. |
| `absent` | TOGA found nothing. **Search the genome directly before calling a loss.** |

Per-mutation details are in `qc.mutations.tsv`, and per-assembly details in
`qc.status_by_assembly.tsv`.

### 2. `qc.master.tsv`: one row per PRDM9 copy
| Column | Meaning |
|---|---|
| `n_zf_array` | Fingers in the C-terminal array |
| `zf0` | 1 if the lone upstream finger (ZF0) is present |
| `KRAB`, `SSXRD`, `SET` | 1 = domain found |
| `synteny`, `synteny_flag` | Share of the usual neighbouring genes found nearby; `NON_SYNTENIC` suggests a copy elsewhere in the genome, not the real gene |
| `reannot_flag` | Flag from step 2 |
| `gate` | `PRDM9` (passes all checks), `PRDM9_like` (has SET but fails something), `not_PRDM9` |
| `gate_reasons` | Why it failed, e.g. `no_SSXRD,non_syntenic` |
| `zf_count_trustworthy` | 1 only for long-read assemblies with a clean ZF exon. **Only use these for ZF counts.** |

### 3. `qc.haplotypes.tsv`: allelic variation and assembly checks
| Column | Meaning |
|---|---|
| `version_check` | Two assembly versions of the **same** haplotype should agree. `DISAGREE` means an assembly problem, not biology. |
| `allelic_zf_diff` | Difference in finger number between hap1 and hap2 of one animal (real within-species variation) |
| `allelic_codes_identical` | 1 if both haplotypes have identical contact residues |

The largest `allelic_zf_diff` is the "noise level": a change on the tree smaller
than this could just be within-species variation.

### 4. `qc.representatives.tsv` and `reps.pass.tsv`
One chosen copy per species. The choice prefers, in order: passes the PRDM9
check, clean re-annotation, best technology, hap1/primary, newest version.
`reps.pass.tsv` keeps only those that pass, and is the set used for all
cross-species analyses.

### Re-annotation report (`prdm9_reannot.report.tsv`, main folder)
| Flag | Meaning |
|---|---|
| `OK` | Clean |
| `GAP_IN_ZF_EXON` | Assembly gap inside the ZF exon: **finger count unreliable** |
| `GAP_NEAR_ZF_EXON` | Gap within 1 kb: be cautious |
| `NO_STOP` | Ran out of sequence before a stop codon |
| `INTERNAL_STOP` | Stop codon inside the protein |
| `NO_START_M` | Protein doesn't start with methionine (start not found) |
| `NO_READTHROUGH_ORF` | Couldn't extend the ZF exon; miniprot's own model used |
| `NO_ALIGNMENT` | Query didn't align to the region |

### ZF array tables
- **`prdm9.per_protein.tsv`:** per copy, giving finger counts, array position,
  `mean_pairwise_id_array` (how similar the fingers are to each other; high
  suggests homogenisation) and `binding_code_string` (contact residues of each
  finger, e.g. `RDER-QSNK-…`).
- **`prdm9.per_zf.tsv`:** one row per finger, with its 28-amino-acid sequence
  and residues at −1, +2, +3, +6. `canonical=0` means unusual spacing.
- **`prdm9.binding_codes.tsv`:** contact residues per finger, ready for clustering.
- **FASTA files:** `prdm9.zf_units.faa/.fna` (each finger separately),
  `prdm9.array.faa` (array only; input for motif predictors),
  `prdm9.nonzf.faa/.fna` (the rest of the protein).

### `prdm9.domains.tsv`
KRAB / SSXRD / SET present (1/0), number of C2H2 hits, `domain_order`, and
`canonical_order` = 1 if the order is KRAB → SSXRD → SET → zinc fingers.

### Selection
- **`zf.repeats.tsv`:** dN/dS among fingers within each array, for `contact`
  vs `noncontact` codons.
  - **Reading the numbers:** dN/dS > 1 suggests positive selection.
  - **The key comparison:** contact vs non-contact within the same arrays.
  - **`frac_boot_contact_gt_noncontact`:** the share of bootstrap replicates
    in which contact > non-contact. Near 1 = consistent.
  - **`note`:** `dS_saturated` means too diverged to interpret; `pS=0` means
    no synonymous changes, so the ratio is undefined.
- **`nonzf.pairs.tsv`:** pairwise dN/dS between species for the non-ZF part.
  `SATURATED` pairs (dS > 1.5) are too far apart to interpret.
- **HyPhy results (`nonzf.busted.json`, `.meme.json`, `.absrel.json`, `.gard.json`):**
  view them by dragging onto the HyPhy Vision website. GARD tests for
  recombination first; if it finds breakpoints, re-run the other tests on the
  partitioned alignment.
  - **BUSTED:** is there selection anywhere in the gene?
  - **MEME:** which sites?
  - **aBSREL:** which branches?

### Concerted evolution
- **`concerted.summary.txt`:** the share of fingers whose closest relative is in
  the same array, compared with what random labelling gives (p-value). High =
  concerted evolution. The "same position in another array" line measures the
  opposite pattern (fingers inherited unchanged by position).
- **`concerted.pairs.tsv`:** for pairs of species in the same family,
  `ratio_within_between` < 1 with small `p_concerted` means fingers are more
  alike within each species than between them. The tests overlap, so report
  them descriptively or correct for multiple testing.
- **`zf_units.treefile`:** tree of all individual fingers; view it in FigTree
  or iTOL and colour by species.

### Tree-based
- **`bats.species.nwk` / `bats.short.nwk`:** the pruned species tree (species
  names / assembly names).
- **`rf_nonzf_vs_species.*`:** Robinson–Foulds distance between the PRDM9 gene
  tree and the species tree.
- **`zfcount.contMap.pdf`:** finger number mapped onto the tree.
- **`zfcount.branch_changes.tsv`:** estimated change in finger number along each
  branch; `exceeds_allelic = TRUE` marks changes bigger than within-species noise.

### Hotspot erosion
- **`erosion.erosion.tsv`:** `erosion` = log2(observed/expected in own genome) −
  median of the same in relatives' genomes. **Negative values** mean the motif is
  rarer in its own genome, as the Red Queen model predicts. `erosion.counts.tsv`
  has the raw counts.

---

## 8. Common problems

| Symptom | Likely cause / fix |
|---|---|
| `command not found: bigBedToBed` | Environment not activated (`mamba activate prdm9`) |
| Many `DOWNLOAD_FAILED` | No internet on compute nodes; run step 1 on a login/data node |
| Many `SEQ_FAILED` | `TWOBIT_PATTERN` is wrong; test with `twoBitInfo` (section 5). After fixing, delete `prdm9_toga_clean/<name>.PRDM9.clean.tsv` for those assemblies and re-run |
| Disk fills up during step 1 | The cache should be deleted per assembly; check `/tmp/udcCache` and `prdm9_tmp/` |
| All statuses look "worst" in QC | TOGA uses different status words; edit `STATUS_RANK` in `prdm9_qc.py` |
| `WARNING: ZF exon unknown` | Set `ZF_EXON` to the number of the ZF exon in the reference transcript |
| `need at least two arrays` | Fewer than two species passed the PRDM9 check; look at `gate_reasons` in `qc.master.tsv` |
| Species missing from the tree | Name differs in the published tree; add to `name_map.tsv` |
| R step fails | Install `ape` and `phytools`; this script has not yet been run on real data |
| `iqtree2: command not found` | See the IQ-TREE note in section 4 |

---

## 9. Limitations

- **Technology labels are guessed.** They come from assembly names; 57 are
  `unknown` until you check them.
- **Finger detection** uses the classic C2H2 pattern (C-x2-C-x12-H-x3-H, allowing
  slightly different spacing). Degenerate fingers (e.g. a missing histidine) are
  not counted. Only canonical fingers enter the dN/dS and concerted-evolution tests.
- **dN/dS among repeats** uses the simple Nei–Gojobori method. Absolute values
  are rough; rely on the **contact vs non-contact comparison**.
- **Motifs:** the erosion test does not predict motifs; you supply them. Long,
  degenerate motifs carry little signal, so use the confident core.
- **One assembly = one allele.** Differences between species can partly reflect
  which allele happened to be sequenced; use `qc.haplotypes.tsv` to judge this.
- **Untested parts:** the pipeline was tested on simulated data. The calls to
  external tools (UCSC tools, miniprot, InterProScan, MAFFT, pal2nal, IQ-TREE,
  HyPhy) and the R script have not been run on real data yet. Check the first
  real run carefully.

---

## Script index

| Script | Step | One-line purpose |
|---|---|---|
| `make_assembly_metadata.py` | 0 | Technology, haplotype and individual for each assembly |
| `run_toga_prdm9.sh` | 1 | Extract PRDM9, neighbours and region DNA; delete the cache |
| `parse_toga_ifne.py` (lab) | 1 | Turn TOGA rows into a readable table |
| `reannotate_prdm9.sh` / `.py` | 2 | Rebuild the protein from DNA with miniprot + read-through ZF exon |
| `prdm9_downstream.sh` | 3 | Runs everything below in order |
| `domain_table.py` | 3.1 | Domain presence and order |
| `prdm9_zf_analysis.py` | 3.2 | Fingers, contact residues, sequence files |
| `prdm9_qc.py` | 3.3 | Status, loss evidence, synteny, PRDM9 check, haplotypes, representatives |
| `fasta_select.py` | 3.4+ | Keep/rename sequences for the representatives |
| `zf_dnds.py` | 3.4, 3.5 | dN/dS among repeats; saturation check |
| `concerted_test.py` | 3.6 | Concerted evolution tests |
| `prune_tree.py` | 3.7 | Prune the published tree to our species |
| `zf_count_asr.R` | 3.7 | Map finger number onto the tree |
| `motif_erosion.py` | 3.8 | Hotspot erosion test |
