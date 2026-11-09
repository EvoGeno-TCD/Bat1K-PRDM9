#!/usr/bin/env python3
"""
dN/dS for PRDM9 with Nei-Gojobori (1986) counting and Jukes-Cantor correction.

Mode 1  repeats  (the Oliver et al. 2009 style test)
    zf_dnds.py repeats <prefix>.zf_units.fna -o zf_dnds
  Within each ZF array, compare every pair of 28-codon repeat units and POOL
  sites and differences over all pairs, separately for
    contact codons     (-1, +2, +3, +6 of the helix)
    non-contact codons (the rest of the unit)
  Pooling is essential: 4 contact codons per pair is far too few on its own.
  CIs come from bootstrapping repeat units within the array.
  Only canonical C-x2-C-x12-H-x3-H units are used so codon positions line up.
  Outputs per array, per species (arrays pooled) and across all arrays.

Mode 2  pairs  (saturation check on a codon alignment, e.g. non-ZF domains)
    zf_dnds.py pairs nonzf.codon.fa -o nonzf_pairs
  Pairwise dN, dS for every pair of sequences; flags dS > --sat (default 1.5)
  or undefined dS (pS >= 0.75), where dN/dS should not be interpreted.

Caveats: NG ignores transition/transversion bias and codon usage, so treat
absolute values as rough; the contact vs non-contact CONTRAST within the same
repeats is the robust quantity. Pairs of repeats are not independent; the
bootstrap over units, not over pairs, is what the CI reflects.
"""
import argparse, itertools, math, random, re, sys
from collections import defaultdict

B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k]
        for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
# contact aa positions within a canonical unit starting at C1:
# C1=0, C2=3, H1=16 -> -1:8  +2:11  +3:12  +6:15
CONTACT_POS = (8, 11, 12, 15)
UNIT = 28
CANON = re.compile(r"^C..C.{12}H...H")

_sites = {}
_diffs = {}


def sites(c):
    """(syn, nonsyn) sites of a codon; stop-changes count as nonsynonymous."""
    if c in _sites:
        return _sites[c]
    s = 0.0
    for p in range(3):
        for b in "ACGT":
            if b != c[p]:
                m = c[:p] + b + c[p + 1:]
                if CODE[m] == CODE[c] and CODE[m] != "*":
                    s += 1 / 3
    _sites[c] = (s, 3 - s)
    return _sites[c]


def diffs(c1, c2):
    """(syn, nonsyn) differences averaged over non-stop pathways."""
    key = (c1, c2)
    if key in _diffs:
        return _diffs[key]
    pos = [i for i in range(3) if c1[i] != c2[i]]
    if not pos:
        return (0.0, 0.0)
    tot_s = tot_n = 0.0
    npath = 0
    for order in itertools.permutations(pos):
        cur, s, n, ok = c1, 0, 0, True
        for p in order:
            nxt = cur[:p] + c2[p] + cur[p + 1:]
            if CODE[nxt] == "*" and nxt != c2:
                ok = False
                break
            if CODE[nxt] == CODE[cur]:
                s += 1
            else:
                n += 1
            cur = nxt
        if ok:
            tot_s += s
            tot_n += n
            npath += 1
    if npath == 0:
        r = (0.0, float(len(pos)))
    else:
        r = (tot_s / npath, tot_n / npath)
    _diffs[key] = r
    return r


def jc(p):
    if p >= 0.75:
        return float("inf")
    return -0.75 * math.log(1 - 4 * p / 3)


def accumulate(seqs, codon_idx):
    """Pool NG counts over all pairs of seqs at the given codon indices."""
    S = N = Sd = Nd = 0.0
    for a, b in itertools.combinations(seqs, 2):
        for ci in codon_idx:
            x, y = a[3 * ci:3 * ci + 3], b[3 * ci:3 * ci + 3]
            if len(x) < 3 or len(y) < 3 or "-" in x + y or "N" in x + y:
                continue
            if CODE.get(x, "*") == "*" or CODE.get(y, "*") == "*":
                continue
            sx, nx = sites(x)
            sy, ny = sites(y)
            S += (sx + sy) / 2
            N += (nx + ny) / 2
            ds, dn = diffs(x, y)
            Sd += ds
            Nd += dn
    return S, N, Sd, Nd


def summarise(S, N, Sd, Nd):
    pS = Sd / S if S else float("nan")
    pN = Nd / N if N else float("nan")
    dS, dN = jc(pS) if S else float("nan"), jc(pN) if N else float("nan")
    w = dN / dS if dS and dS not in (float("inf"),) and not math.isnan(dS) else float("nan")
    return pN, pS, dN, dS, w


def fmt(x):
    if isinstance(x, float) and x == 0:
        x = 0.0
    return "NA" if x is None or (isinstance(x, float) and (math.isnan(x))) else (
        "inf" if x == float("inf") else f"{x:.4f}")


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line.upper())
    return {k: "".join(v) for k, v in out.items()}


def translate(nt):
    return "".join(CODE.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def repeats_mode(a):
    units = read_fasta(a.fasta)
    arrays = defaultdict(list)
    skipped = 0
    for name, s in units.items():
        if len(s) < 3 * UNIT or not CANON.match(translate(s[:3 * UNIT])):
            skipped += 1
            continue
        arrays[name.rsplit("|ZF", 1)[0]].append(s[:3 * UNIT])
    print(f"{sum(map(len, arrays.values()))} canonical units in {len(arrays)} arrays "
          f"({skipped} non-canonical skipped)", file=sys.stderr)
    contact = list(CONTACT_POS)
    noncontact = [i for i in range(UNIT) if i not in CONTACT_POS]
    rng = random.Random(a.seed)

    def analyse(seqs):
        res = {}
        for cls, idx in (("contact", contact), ("noncontact", noncontact)):
            res[cls] = summarise(*accumulate(seqs, idx))
        # bootstrap ratio (contact w / noncontact w) over units
        boots = []
        for _ in range(a.boot):
            bs = [rng.choice(seqs) for _ in seqs]
            wc = summarise(*accumulate(bs, contact))[4]
            wn = summarise(*accumulate(bs, noncontact))[4]
            boots.append((wc, wn))
        return res, boots

    def ci(vals):
        v = sorted(x for x in vals if not math.isnan(x) and x != float("inf"))
        if len(v) < 10:
            return ("NA", "NA")
        return fmt(v[int(0.025 * len(v))]), fmt(v[int(0.975 * len(v)) - 1])

    hdr = ("level\tname\tn_units\t" + "\t".join(
        f"{c}_{k}" for c in ("contact", "noncontact") for k in ("pN", "pS", "dN", "dS", "dNdS"))
        + "\tcontact_dNdS_CI\tnoncontact_dNdS_CI\tfrac_boot_contact_gt_noncontact\tnote\n")
    with open(a.out + ".repeats.tsv", "w") as fh:
        fh.write(hdr)
        groups = [("array", k, v) for k, v in sorted(arrays.items())]
        by_sp = defaultdict(list)
        for k, v in arrays.items():
            by_sp[k.split("|")[0]].append(v)
        # species level: pool pairs within each array, never across arrays
        for level, name, seqs in groups:
            if len(seqs) < 2:
                continue
            res, boots = analyse(seqs)
            write(fh, level, name, len(seqs), res, boots, ci)
        for sp, arrs in sorted(by_sp.items()):
            write_pooled(fh, "species", sp, arrs, contact, noncontact)
        write_pooled(fh, "all", "all_arrays", list(arrays.values()), contact, noncontact)


def write(fh, level, name, n, res, boots, ci):
    c, nc = res["contact"], res["noncontact"]
    wc = [b[0] for b in boots]
    wn = [b[1] for b in boots]
    valid = [(x, y) for x, y in boots if not (math.isnan(x) or math.isnan(y))]
    frac = fmt(sum(x > y for x, y in valid) / len(valid)) if valid else "NA"
    note = []
    if c[3] == float("inf") or nc[3] == float("inf") or (not math.isnan(nc[3]) and nc[3] > 1.5):
        note.append("dS_saturated")
    if (not math.isnan(c[1]) and c[1] == 0) or (not math.isnan(nc[1]) and nc[1] == 0):
        note.append("pS=0")
    fh.write(f"{level}\t{name}\t{n}\t" + "\t".join(fmt(x) for x in c + nc)
             + f"\t{'-'.join(ci(wc))}\t{'-'.join(ci(wn))}\t{frac}\t{','.join(note) or 'ok'}\n")


def write_pooled(fh, level, name, arrs, contact, noncontact):
    tot = {}
    for cls, idx in (("contact", contact), ("noncontact", noncontact)):
        acc = [0.0] * 4
        for seqs in arrs:
            if len(seqs) >= 2:
                acc = [x + y for x, y in zip(acc, accumulate(seqs, idx))]
        tot[cls] = summarise(*acc)
    n = sum(len(s) for s in arrs)
    c, nc = tot["contact"], tot["noncontact"]
    fh.write(f"{level}\t{name}\t{n}\t" + "\t".join(fmt(x) for x in c + nc)
             + "\tNA\tNA\tNA\tpooled_within_arrays\n")


def pairs_mode(a):
    seqs = read_fasta(a.fasta)
    names = list(seqs)
    L = min(len(s) for s in seqs.values()) // 3
    with open(a.out + ".pairs.tsv", "w") as fh:
        fh.write("seq1\tseq2\tcodons_compared\tpN\tpS\tdN\tdS\tdNdS\tflag\n")
        for x, y in itertools.combinations(names, 2):
            sx, sy = seqs[x], seqs[y]
            idx = [i for i in range(L) if "-" not in sx[3*i:3*i+3] + sy[3*i:3*i+3]]
            pN, pS, dN, dS, w = summarise(*accumulate([sx, sy], idx))
            flag = "SATURATED" if (dS == float("inf") or (not math.isnan(dS) and dS > a.sat)) else "ok"
            fh.write(f"{x}\t{y}\t{len(idx)}\t{fmt(pN)}\t{fmt(pS)}\t{fmt(dN)}\t{fmt(dS)}\t{fmt(w)}\t{flag}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["repeats", "pairs"])
    ap.add_argument("fasta")
    ap.add_argument("-o", "--out", default="zf_dnds")
    ap.add_argument("--boot", type=int, default=200)
    ap.add_argument("--sat", type=float, default=1.5)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    repeats_mode(a) if a.mode == "repeats" else pairs_mode(a)


if __name__ == "__main__":
    main()
