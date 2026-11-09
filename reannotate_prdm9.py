#!/usr/bin/env python3
"""
Rebuild PRDM9 CDS/protein from the genome, with the ZF exon read directly
from the sequence instead of trusting a projection/alignment.

For each locus FASTA (header  >species|chrom|offset  as written by
run_toga_prdm9.sh) and its miniprot GFF:
  1. take miniprot's best model (Rank=1, highest score);
  2. find the most upstream CDS segment k from which an uninterrupted
     open reading frame (same frame as the model) runs through ALL later
     segments - this repairs fake introns that aligners put inside the
     repetitive ZF exon;
  3. read on from segment k to the first in-frame stop (so arrays longer
     than the query are captured in full), capped at --max-exon bp;
  4. report N/gap content in that exon - gaps inside or next to the array
     mean the ZF count is not trustworthy.

usage: reannotate_prdm9.py --gff-dir gff/ --seq-dir prdm9_locus_seq/ -o reannot
outputs: reannot.faa reannot.fna reannot.report.tsv
"""
import argparse, glob, os, re, sys

CODON = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            CODON[a + b + c] = _aa[16 * i + 4 * j + k]
ZF_RE = re.compile(r"C.{2,4}C.{12}H.{3,5}H")
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def translate(nt):
    return "".join(CODON.get(nt[i:i + 3].upper(), "X") for i in range(0, len(nt) - 2, 3))


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seq.append(line.strip())
    return name, "".join(seq).upper()


def parse_gff(path):
    """-> {mrna_id: dict(seqid, strand, score, rank, identity, target, cds=[(s,e,phase)])}"""
    models = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            models[attr["ID"]] = dict(
                seqid=f[0], strand=f[6], score=float(f[5]) if f[5] != "." else 0.0,
                rank=int(attr.get("Rank", 1)), identity=attr.get("Identity", "NA"),
                target=attr.get("Target", "NA").split()[0], cds=[])
        elif f[2] == "CDS" and attr.get("Parent") in models:
            models[attr["Parent"]]["cds"].append(
                (int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] in "012" else 0))
    return models


def rebuild(genome, m, max_exon):
    """Return dict with nt, aa and QC for one miniprot model."""
    L = len(genome)
    if m["strand"] == "-":
        seq = rc(genome)
        segs = sorted(((L - e, L - s, p) for s, e, p in m["cds"]), key=lambda x: x[0])
    else:
        seq = genome
        segs = sorted(m["cds"], key=lambda x: x[0])
    phase0 = segs[0][2]
    parts = [seq[s:e] for s, e, _ in segs]
    starts_t, t = [], -phase0            # transcript coord (after trimming phase0)
    for s, e, _ in segs:
        starts_t.append(t)
        t += e - s
    model_nt = "".join(parts)[phase0:]

    def readthrough(k):
        Tk = max(starts_t[k], 0)
        gk = segs[k][0] + (Tk - starts_t[k])
        prefix = model_nt[:Tk]
        cand = prefix + seq[gk:gk + max_exon + 3]
        stop_at = None
        for ci in range(0, len(cand) - 2, 3):
            if CODON.get(cand[ci:ci + 3], "X") == "*" and ci + 3 > Tk:
                stop_at = ci + 3
                break
        end_g = gk + ((stop_at if stop_at else len(cand)) - Tk)
        return (k, prefix, gk, stop_at, end_g), end_g >= segs[-1][1] - 2

    # Walk upstream from the last segment. Merge segment k into the read-
    # through exon only if (a) the ORF from k runs through everything
    # downstream and (b) the gap after k looks like array, not intron:
    # it is tiny (<= 90 bp) or its in-frame translation contains a C2H2.
    chosen, _ = readthrough(len(segs) - 1)
    for k in range(len(segs) - 2, -1, -1):
        cand, ok = readthrough(k)
        if not ok:
            break
        gap_s, gap_e = segs[k][1], segs[k + 1][0]
        gap_nt = seq[gap_s:gap_e]
        Tk = max(starts_t[k], 0)
        off = (gap_s - cand[2] + Tk) % 3          # transcript offset of gap start
        gap_aa = translate(gap_nt[(3 - off) % 3:])
        if (gap_e - gap_s) <= 90 or ZF_RE.search(gap_aa):
            chosen = cand
        else:
            break
    if not chosen[3] and chosen[4] < segs[-1][1] - 2:
        chosen = None
    if chosen is None:                        # fall back to miniprot's own model
        aa = translate(model_nt)
        return dict(nt=model_nt, aa=aa.split("*")[0], k=None, stop=aa.find("*") >= 0,
                    zf_exon=(None, None), n_in_exon=None, changed=0, ok=0)
    k, prefix, gk, stop_at, end_g = chosen
    nt = (prefix + seq[gk:end_g]) if stop_at else (prefix + seq[gk:gk + max_exon])
    nt = nt[:len(nt) - len(nt) % 3]
    aa = translate(nt).rstrip("*")
    exon = seq[gk:end_g]
    flank = seq[max(gk - 1000, 0):min(end_g + 1000, L)]
    # back to forward-strand local coords for reporting
    zs, ze = (gk, end_g) if m["strand"] == "+" else (L - end_g, L - gk)
    return dict(nt=nt, aa=aa, k=k, stop=bool(stop_at), zf_exon=(zs, ze),
                n_in_exon=exon.count("N"), n_flank=flank.count("N"),
                changed=int(k < len(segs) - 1 or end_g > segs[-1][1]), ok=1,
                n_segs=len(segs))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--seq-dir", required=True)
    ap.add_argument("-o", "--out", default="reannot")
    ap.add_argument("--max-exon", type=int, default=6000,
                    help="max length of the read-through ZF exon (bp)")
    args = ap.parse_args()

    faa, fna = open(args.out + ".faa", "w"), open(args.out + ".fna", "w")
    rep = open(args.out + ".report.tsv", "w")
    rep.write("id\tspecies\tlocus_fasta\tquery\tidentity\tstrand\tn_cds_segments\t"
              "zf_exon_seg\tzf_exon_start\tzf_exon_end\tzf_exon_len\tstop_found\t"
              "N_in_zf_exon\tN_within_1kb\tmodel_changed\tprot_len\tstarts_with_M\tflag\n")
    for fa in sorted(glob.glob(os.path.join(args.seq_dir, "*.fa"))):
        base = os.path.basename(fa)[:-3]
        gff = os.path.join(args.gff_dir, base + ".gff")
        if not os.path.exists(gff):
            print(f"no GFF for {base}", file=sys.stderr)
            continue
        name, genome = read_fasta(fa)
        sp, chrom, off = name.split("|")
        off = int(off)
        models = [m for m in parse_gff(gff).values() if m["cds"]]
        if not models:
            rep.write(f"NA\t{sp}\t{base}\t" + "NA\t" * 14 + "NO_ALIGNMENT\n")
            continue
        m = max(models, key=lambda x: (x["rank"] == 1, x["score"]))
        r = rebuild(genome, m, args.max_exon)
        zs, ze = r["zf_exon"]
        gs, ge = (off + zs, off + ze) if zs is not None else ("NA", "NA")
        pid = f"{sp}|{chrom}:{gs}-{ge}"
        flags = []
        if not r["ok"]:
            flags.append("NO_READTHROUGH_ORF")
        if not r["stop"]:
            flags.append("NO_STOP")
        if r.get("n_in_exon"):
            flags.append("GAP_IN_ZF_EXON")
        elif r.get("n_flank"):
            flags.append("GAP_NEAR_ZF_EXON")
        if "*" in r["aa"]:
            flags.append("INTERNAL_STOP")
        if not r["aa"].startswith("M"):
            flags.append("NO_START_M")
        faa.write(f">{pid}\n{r['aa']}\n")
        fna.write(f">{pid}\n{r['nt']}\n")
        rep.write("\t".join(map(str, [
            pid, sp, base, m["target"], m["identity"], m["strand"], len(m["cds"]),
            r["k"] if r["k"] is not None else "NA", gs, ge,
            (ge - gs) if zs is not None else "NA", int(r["stop"]),
            r.get("n_in_exon", "NA"), r.get("n_flank", "NA"), r["changed"],
            len(r["aa"]), int(r["aa"].startswith("M")), ",".join(flags) or "OK"])) + "\n")
    for f in (faa, fna, rep):
        f.close()


if __name__ == "__main__":
    main()
