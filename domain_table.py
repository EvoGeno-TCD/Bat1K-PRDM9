#!/usr/bin/env python3
"""
Domain presence per protein from InterProScan TSV.
  domain_table.py prdm9.interpro.tsv prdm9_reannot.faa > prdm9.domains.tsv

Calls by signature accession first, falling back to whole-word matches in the
descriptions. Also reports domain ORDER (PRDM9 = KRAB -> SSXRD -> SET -> ZFs)
and proteins with no hits at all (listed with zeros, so losses are visible).
"""
import re, sys
from collections import defaultdict

ACC = {
    "KRAB": {"PF01352", "SM00349", "PS50805"},
    "SSXRD": {"PF09514"},
    "SET": {"PF00856", "SM00317", "PS50280"},
    "C2H2": {"PF00096", "PF13465", "PF13912", "SM00355", "PS50157", "PS00028"},
}
TXT = {
    "KRAB": re.compile(r"\bkrab\b|kruppel-associated box", re.I),
    "SSXRD": re.compile(r"\bssxrd\b|ssx repression", re.I),
    "SET": re.compile(r"\bset domain\b|\bpr domain\b|\bpr[-/]set\b|\bpr/set\b", re.I),
    "C2H2": re.compile(r"\bc2h2\b|zinc finger,? c2h2", re.I),
}


def main(ipr, faa):
    ids = [l[1:].split()[0] for l in open(faa) if l.startswith(">")]
    hits = defaultdict(lambda: defaultdict(list))
    for line in open(ipr):
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        pid, acc = f[0], f[4]
        desc = " ".join(f[5:6] + (f[12:13] if len(f) > 12 else []))
        start = int(f[6])
        for dom in ACC:
            if acc in ACC[dom] or TXT[dom].search(desc):
                hits[pid][dom].append(start)
    print("id\tKRAB\tSSXRD\tSET\tn_C2H2_hits\tdomain_order\tcanonical_order")
    for pid in ids:
        h = hits.get(pid, {})
        first = {d: min(v) for d, v in h.items() if v}
        order = "-".join(d for d, _ in sorted(first.items(), key=lambda x: x[1]))
        canon = int(order.startswith("KRAB-SSXRD-SET") and "C2H2" in order)
        n_c2h2 = len(set(h.get("C2H2", [])))
        print(f"{pid}\t{int('KRAB' in first)}\t{int('SSXRD' in first)}\t{int('SET' in first)}\t"
              f"{n_c2h2}\t{order or 'none'}\t{canon}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
