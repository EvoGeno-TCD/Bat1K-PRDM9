#!/usr/bin/env python3
"""
Hotspot-erosion ('Red Queen') test: is each species' predicted PRDM9 motif
depleted in its OWN genome relative to the genomes of other species?

Inputs
  motifs.tsv   species<TAB>motif   (IUPAC consensus, e.g. CCNCCNTNNCCNC)
               - one motif per species/allele; get these from a C2H2
               recognition-code predictor run on <prefix>.array.faa
               (e.g. the Persikov & Singh predictor). This script does NOT
               predict motifs.
  genomes.tsv  species<TAB>/path/genome.fa[.gz]   (same species labels)

For every motif x genome: observed count (both strands, overlapping) and
expected count under an order-1 Markov model of that genome (keeps GC and
CpG depletion, which matter a lot for GC-rich PRDM9 motifs).
  oe = observed / expected
  erosion = log2(oe in own genome) - median log2(oe in the other genomes)
Negative erosion = motif is rarer in its own genome than elsewhere, the
signature expected if PRDM9-directed gene conversion has been destroying it.

Caveats: long motifs with many N's have tiny information; use the core
(highest-confidence) part of the predicted motif. Compare close relatives
(e.g. --meta + --same-family) - distant genomes differ in composition more
than any Markov-1 model fixes. Masked (lowercase) sequence is counted.
"""
import argparse, csv, gzip, math, re, statistics, sys
import numpy as np

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG",
         "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT",
         "V": "ACG", "N": "ACGT"}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W",
        "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}
IDX = {"A": 0, "C": 1, "G": 2, "T": 3}
LUT = np.full(256, 4, dtype=np.uint8)
for b, i in IDX.items():
    LUT[ord(b)] = i
    LUT[ord(b.lower())] = i


def rc_motif(m):
    return "".join(COMP[c] for c in reversed(m))


def to_regex(m):
    return re.compile("(?=" + "".join(f"[{IUPAC[c]}]" if len(IUPAC[c]) > 1 else c for c in m) + ")")


def chroms(path):
    op = gzip.open if path.endswith(".gz") else open
    name, buf = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if buf:
                    yield "".join(buf)
                buf = []
            else:
                buf.append(line.strip())
    if buf:
        yield "".join(buf)


def expected(m, p1, P):
    """Expected count per valid window under Markov-1, via DP over states."""
    v = np.array([p1[i] if "ACGT"[i] in IUPAC[m[0]] else 0 for i in range(4)])
    for c in m[1:]:
        mask = np.array([1.0 if b in IUPAC[c] else 0.0 for b in "ACGT"])
        v = (v @ P) * mask
    return v.sum()


def scan_genome(path, motifs):
    di = np.zeros((5, 5))
    obs = {m: 0 for m in motifs}
    pats = {m: [to_regex(m)] + ([to_regex(rc_motif(m))] if rc_motif(m) != m else [])
            for m in motifs}
    total = 0
    for seq in chroms(path):
        codes = LUT[np.frombuffer(seq.encode(), dtype=np.uint8)]
        _fast_di(di, codes)
        up = seq.upper()
        total += int((codes < 4).sum())
        for m, ps in pats.items():
            for p in ps:
                obs[m] += sum(1 for _ in p.finditer(up))
    D = di[:4, :4]
    P = D / D.sum(axis=1, keepdims=True)
    p1 = D.sum(axis=1) / D.sum()
    out = {}
    for m in motifs:
        e = total * (expected(m, p1, P) + (expected(rc_motif(m), p1, P)
                                          if rc_motif(m) != m else 0))
        out[m] = (obs[m], e)
    return out, total


def _fast_di(di, codes):
    pair = codes[:-1].astype(np.int64) * 5 + codes[1:]
    di += np.bincount(pair, minlength=25).reshape(5, 5)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("motifs")
    ap.add_argument("genomes")
    ap.add_argument("--meta", help="bat_assemblies.tsv (species -> family)")
    ap.add_argument("--same-family", action="store_true",
                    help="compare each motif only against genomes of the same family")
    ap.add_argument("-o", "--out", default="erosion")
    a = ap.parse_args()

    motifs = [l.rstrip("\n").split("\t")[:2] for l in open(a.motifs) if l.strip() and not l.startswith("#")]
    genomes = [l.rstrip("\n").split("\t")[:2] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    fam = {}
    if a.meta:
        for r in csv.DictReader(open(a.meta), delimiter="\t"):
            fam[r["species"]] = r["family"]
            fam[r["short"]] = r["family"]
    uniq = sorted({m.upper() for _, m in motifs})
    oe = {}
    with open(a.out + ".counts.tsv", "w") as fh:
        fh.write("genome_species\tmotif\tobserved\texpected\tobs_exp\tvalid_bases\n")
        for gsp, path in genomes:
            need = uniq if not a.same_family else sorted(
                {m.upper() for sp, m in motifs if fam.get(sp) == fam.get(gsp)})
            if not need:
                continue
            print(f"scanning {gsp} for {len(need)} motifs", file=sys.stderr)
            res, tot = scan_genome(path, need)
            for m, (o, e) in res.items():
                r = o / e if e > 0 else float("nan")
                oe[(gsp, m)] = r
                fh.write(f"{gsp}\t{m}\t{o}\t{e:.1f}\t{r:.4f}\t{tot}\n")
    with open(a.out + ".erosion.tsv", "w") as fh:
        fh.write("motif_species\tmotif\tlog2oe_own\tmedian_log2oe_others\tn_others\terosion\t"
                 "rank_own_among_all\n")
        gsp_all = [g for g, _ in genomes]
        for sp, m in motifs:
            m = m.upper()
            if (sp, m) not in oe:
                fh.write(f"{sp}\t{m}\tNA\tNA\t0\tNA\tNA\n")
                continue
            others = [g for g in gsp_all if g != sp and
                      (not a.same_family or fam.get(g) == fam.get(sp))]
            lo = [math.log2(oe[(g, m)]) for g in others if (oe.get((g, m)) or 0) > 0]
            own = math.log2(oe[(sp, m)]) if oe[(sp, m)] > 0 else float("-inf")
            med = statistics.median(lo) if lo else float("nan")
            allv = sorted([own] + lo)
            fh.write(f"{sp}\t{m}\t{own:.3f}\t{med:.3f}\t{len(lo)}\t{own - med:.3f}\t"
                     f"{allv.index(own) + 1}/{len(allv)}\n")


if __name__ == "__main__":
    main()
