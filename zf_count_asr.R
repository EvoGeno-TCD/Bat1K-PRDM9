#!/usr/bin/env Rscript
# Map PRDM9 ZF-array length onto the bat tree and flag lineage-specific
# expansions/contractions that exceed within-individual allelic variation.
# NOTE: written without access to R for testing - check the first run.
#
# Rscript zf_count_asr.R bats.short.nwk prdm9_qc.representatives.tsv \
#         prdm9_qc.master.tsv prdm9_qc.haplotypes.tsv out_prefix
suppressPackageStartupMessages({ library(ape); library(phytools) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: tree reps master haplotypes out_prefix")
tree <- read.tree(args[1])
reps <- read.delim(args[2], stringsAsFactors = FALSE)
master <- read.delim(args[3], stringsAsFactors = FALSE)
haps <- read.delim(args[4], stringsAsFactors = FALSE)
out <- args[5]

# keep representatives that pass the PRDM9 gate and have a trustworthy count
ok <- master$id[master$gate == "PRDM9" & master$zf_count_trustworthy == 1]
reps <- reps[reps$id %in% ok, ]
x <- setNames(as.numeric(reps$n_zf_array), reps$short)
x <- x[!is.na(x)]
tree <- drop.tip(tree, setdiff(tree$tip.label, names(x)))
x <- x[tree$tip.label]
cat(sprintf("%d species with trustworthy ZF counts\n", length(x)))

# allelic noise: the largest hap1-vs-hap2 difference seen in any individual
adiff <- suppressWarnings(as.numeric(haps$allelic_zf_diff))
thr <- if (all(is.na(adiff))) 1 else max(adiff, na.rm = TRUE)
cat(sprintf("allelic threshold (max hap1/hap2 ZF difference): %g\n", thr))

# phylogenetic signal
K <- phylosig(tree, x, method = "K", test = TRUE)
L <- phylosig(tree, x, method = "lambda", test = TRUE)
cat(sprintf("Blomberg K = %.3f (p = %.3g); Pagel lambda = %.3f (p = %.3g)\n",
            K$K, K$P, L$lambda, L$P))

# ancestral states (Brownian motion on counts - a reasonable first pass;
# for small integer counts also consider an ordered Mk model / fitDiscrete)
anc <- fastAnc(tree, x, CI = TRUE)
vals <- c(x, anc$ace)            # tips 1..n, then internal nodes n+1..
names(vals)[seq_along(x)] <- seq_along(x)
e <- tree$edge
delta <- vals[as.character(e[, 2])] - vals[as.character(e[, 1])]
desc <- sapply(e[, 2], function(n) {
  if (n <= Ntip(tree)) tree$tip.label[n]
  else paste(extract.clade(tree, n)$tip.label, collapse = ",")
})
br <- data.frame(parent = e[, 1], child = e[, 2],
                 parent_state = round(vals[as.character(e[, 1])], 2),
                 child_state = round(vals[as.character(e[, 2])], 2),
                 change = round(delta, 2),
                 exceeds_allelic = abs(delta) > thr,
                 descendants = desc)
br <- br[order(-abs(br$change)), ]
write.table(br, paste0(out, ".branch_changes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(paste0(out, ".contMap.pdf"), width = 8, height = max(6, Ntip(tree) * 0.12))
cm <- contMap(tree, x, plot = FALSE)
plot(setMap(cm, c("#2166ac", "#f7f7f7", "#b2182b")), fsize = c(0.5, 0.8),
     leg.txt = "ZFs in array", lwd = 3)
dev.off()
cat("wrote", paste0(out, ".branch_changes.tsv"), "and", paste0(out, ".contMap.pdf"), "\n")
