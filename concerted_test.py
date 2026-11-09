#!/usr/bin/env python3
"""
Test for concerted evolution of PRDM9 ZF repeats.

Input: <prefix>.zf_units.faa (or .fna) from prdm9_zf_analysis.py, IDs like
       <assembly>|<locus>|ZF07.  Restrict to one array per species with
       --keep representatives.tsv (column 'id'), otherwise all arrays are used.

Two complementary analyses:

1. Nearest neighbours. For each repeat, find its closest other repeat
   (p-distance; ties split evenly) and record whether it lies in
     - the SAME array                (concerted evolution / recent duplication)
     - the SAME POSITION of another species' array (positional orthology)
   Null: species labels permuted among repeats, keeping how many repeats each
   array has (--perm times). Reports observed fraction, null mean and p.

2. Species pairs. For each pair of arrays A, B (within --pairs file, or all
   pairs in the same family via --meta, or all pairs):
     within = mean distance among repeats inside A and inside B
     between = mean distance between repeats of A and repeats of B
   Concerted evolution predicts within < between even for close relatives
   (repeats homogenised within each lineage). Positional orthology with old
   repeats predicts within >= between. p-value: permute the A/B labels across
   the pooled repeats (sizes kept) and ask how often within/between is as
   small as observed.

Output: <o>.nn.tsv, <o>.pairs.tsv, <o>.summary.txt, <o>.dist.tsv (optional matrix)
"""
import argparse, csv, itertools, random, sys
from collections import defaultdict
import numpy as np


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


def pdist_matrix(seqs):
    L = min(len(s) for s in seqs)
    arr = np.array([np.frombuffer(s[:L].encode(), dtype=np.uint8) for s in seqs])
    n = len(seqs)
    D = np.zeros((n, n))
    for i in range(n):
        D[i] = (arr != arr[i]).mean(axis=1)
    return D


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("units")
    ap.add_argument("--keep", help="representatives.tsv (column id) to use one array per species")
    ap.add_argument("--meta", help="bat_assemblies.tsv: restrict pair tests to same family")
    ap.add_argument("--pairs", help="TSV of assembly pairs to test (two columns)")
    ap.add_argument("--perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--write-matrix", action="store_true")
    ap.add_argument("-o", "--out", default="concerted")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    units = read_fasta(a.units)
    keep = None
    if a.keep:
        keep = {r["id"] for r in csv.DictReader(open(a.keep), delimiter="\t")}
    names, seqs, arr_of, pos_of = [], [], [], []
    for n, s in units.items():
        arr, zf = n.rsplit("|ZF", 1)
        if keep is not None and arr not in keep:
            continue
        names.append(n); seqs.append(s); arr_of.append(arr); pos_of.append(int(zf))
    arrays = sorted(set(arr_of))
    if len(arrays) < 2:
        sys.exit("need at least two arrays")
    aidx = np.array([arrays.index(x) for x in arr_of])
    pos = np.array(pos_of)
    D = pdist_matrix(seqs)
    np.fill_diagonal(D, np.inf)
    print(f"{len(names)} repeats from {len(arrays)} arrays", file=sys.stderr)

    # ---- 1. nearest neighbours ------------------------------------------
    def nn_stats(lab):
        same_arr = same_pos = 0.0
        for i in range(len(lab)):
            m = D[i].min()
            nb = np.where(D[i] == m)[0]
            same_arr += np.mean(lab[nb] == lab[i])
            same_pos += np.mean((lab[nb] != lab[i]) & (pos[nb] == pos[i]))
        return same_arr / len(lab), same_pos / len(lab)

    obs_arr, obs_pos = nn_stats(aidx)
    null_arr = []
    for _ in range(a.perm):
        null_arr.append(nn_stats(rng.permutation(aidx))[0])
    null_arr = np.array(null_arr)
    p_arr = (1 + np.sum(null_arr >= obs_arr)) / (1 + a.perm)

    with open(a.out + ".nn.tsv", "w") as fh:
        fh.write("repeat\tarray\tnearest\tdistance\tnearest_same_array\n")
        for i, n in enumerate(names):
            j = int(np.argmin(D[i]))
            fh.write(f"{n}\t{arr_of[i]}\t{names[j]}\t{D[i, j]:.3f}\t{int(aidx[j] == aidx[i])}\n")

    # ---- 2. species pairs ------------------------------------------------
    fam = {}
    if a.meta:
        for r in csv.DictReader(open(a.meta), delimiter="\t"):
            fam[r["short"]] = r["family"]
    asm_of = {x: x.split("|")[0] for x in arrays}
    if a.pairs:
        wanted = {tuple(sorted(l.split()[:2])) for l in open(a.pairs) if l.strip()}
        tests = [(x, y) for x, y in itertools.combinations(arrays, 2)
                 if tuple(sorted((asm_of[x], asm_of[y]))) in wanted]
    elif fam:
        tests = [(x, y) for x, y in itertools.combinations(arrays, 2)
                 if fam.get(asm_of[x]) and fam.get(asm_of[x]) == fam.get(asm_of[y])]
    else:
        tests = list(itertools.combinations(arrays, 2))

    Dm = D.copy()
    np.fill_diagonal(Dm, np.nan)

    def wb(ix, lab):
        A, Bm = ix[lab], ix[~lab]
        w = np.concatenate([Dm[np.ix_(A, A)][np.triu_indices(len(A), 1)],
                            Dm[np.ix_(Bm, Bm)][np.triu_indices(len(Bm), 1)]])
        b = Dm[np.ix_(A, Bm)].ravel()
        return np.nanmean(w), np.nanmean(b)

    n_sig = 0
    with open(a.out + ".pairs.tsv", "w") as fh:
        fh.write("array_A\tarray_B\tnA\tnB\twithin\tbetween\tratio_within_between\tp_concerted\n")
        for x, y in tests:
            ia = np.where(aidx == arrays.index(x))[0]
            ib = np.where(aidx == arrays.index(y))[0]
            if len(ia) < 2 or len(ib) < 2:
                continue
            ix = np.concatenate([ia, ib])
            lab = np.array([True] * len(ia) + [False] * len(ib))
            w, b = wb(ix, lab)
            r = w / b if b > 0 else np.nan
            null = []
            for _ in range(a.perm):
                pw, pb = wb(ix, rng.permutation(lab))
                null.append(pw / pb if pb > 0 else np.nan)
            null = np.array(null)
            p = (1 + np.nansum(null <= r)) / (1 + np.sum(~np.isnan(null)))
            n_sig += p < 0.05
            fh.write(f"{x}\t{y}\t{len(ia)}\t{len(ib)}\t{w:.3f}\t{b:.3f}\t{r:.3f}\t{p:.4f}\n")

    with open(a.out + ".summary.txt", "w") as fh:
        fh.write(f"repeats: {len(names)}  arrays: {len(arrays)}\n")
        fh.write(f"nearest neighbour in same array: {obs_arr:.3f} "
                 f"(null mean {null_arr.mean():.3f}, p = {p_arr:.4f})\n")
        fh.write(f"nearest neighbour at same position in another array (orthology signal): {obs_pos:.3f}\n")
        fh.write(f"species-pair tests: {len(tests)}; within<between at p<0.05: {n_sig}\n")
        fh.write("Note: pair tests are not independent; report them descriptively or "
                 "correct for multiple testing.\n")
    print(open(a.out + ".summary.txt").read(), file=sys.stderr)

    if a.write_matrix:
        with open(a.out + ".dist.tsv", "w") as fh:
            fh.write("\t" + "\t".join(names) + "\n")
            Dw = D.copy(); np.fill_diagonal(Dw, 0)
            for i, n in enumerate(names):
                fh.write(n + "\t" + "\t".join(f"{v:.3f}" for v in Dw[i]) + "\n")


if __name__ == "__main__":
    main()
