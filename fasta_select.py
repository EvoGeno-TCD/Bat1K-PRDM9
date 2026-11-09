#!/usr/bin/env python3
"""
Subset (and optionally rename) FASTA records using a representatives table.
  fasta_select.py in.fa representatives.tsv [--rename] [--units] > out.fa
--rename  rename records to the 'short' column (Newick/HyPhy-safe names)
--units   records are ZF units (<id>|ZFnn): match on <id>, keep the |ZFnn suffix
"""
import csv, sys

args = [a for a in sys.argv[1:] if not a.startswith("--")]
rename, units = "--rename" in sys.argv, "--units" in sys.argv
if len(args) != 2:
    sys.exit(__doc__)
rep = {r["id"]: r["short"] for r in csv.DictReader(open(args[1]), delimiter="\t")}
keep = False
for line in open(args[0]):
    if line.startswith(">"):
        name = line[1:].split()[0]
        base, suf = (name.rsplit("|ZF", 1) + [""])[:2] if units else (name, "")
        keep = base in rep
        if keep:
            new = rep[base] if rename else base
            sys.stdout.write(f">{new}{'_ZF' + suf if (units and rename) else ('|ZF' + suf if units else '')}\n")
    elif keep:
        sys.stdout.write(line)
