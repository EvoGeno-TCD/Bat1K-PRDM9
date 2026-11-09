#!/usr/bin/env python3
"""
PRDM9 zinc-finger array dissection.

Input : protein FASTA of PRDM9 candidates (one record per copy/haplotype),
        optionally the matching CDS FASTA (same IDs) for codon-level outputs.
        IDs should look like  <species>|<copy_id>  (anything before the first
        '|' is used as the species/assembly label).

Output (prefix -o):
  <o>.per_protein.tsv   ZF counts, array length, ZF0 presence, array identity
  <o>.per_zf.tsv        every ZF: position, 28-aa unit, contact residues
  <o>.zf_units.faa      one record per array ZF unit (for repeat trees)
  <o>.zf_units.fna      same, nucleotide (if --cds)
  <o>.zf_contact.fna    contact codons only per unit  (if --cds)
  <o>.zf_noncontact.fna non-contact codons per unit    (if --cds)
  <o>.nonzf.faa         protein with the terminal ZF array removed (for
                        aligning KRAB/SSXRD/SET across species)
  <o>.nonzf.fna         matching CDS (if --cds) -> pal2nal / HyPhy
  <o>.array.faa         the ZF array only (input for ZF binding predictors)
  <o>.binding_codes.tsv contact-residue string per ZF, for clustering

C2H2 helix numbering: with the first zinc-coordinating His at helix
position +7, the DNA-contact residues -1, +2, +3, +6 sit at H1-8, H1-5,
H1-4 and H1-1 (Zif268 finger 1 'RSDELTRH' -> R D E R).
"""
import argparse, itertools, re, sys
from collections import defaultdict

ZF_RE = re.compile(r"(?=(C.{2,4}C.{12}H.{3,5}H))")   # overlapping scan
CONTACT = {"-1": -8, "+2": -5, "+3": -4, "+6": -1}    # offsets from H1
UNIT_LEN = 28
ARRAY_MAX_GAP = 34          # max C1-to-C1 distance for tandem neighbours


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {k: "".join(v).upper().rstrip("*") for k, v in seqs.items()}


def find_zfs(prot):
    hits = []
    for m in ZF_RE.finditer(prot):
        motif = m.group(1)
        c1 = m.start()
        mm = re.match(r"C(.{2,4})C(.{12})H(.{3,5})H", motif)
        c2 = c1 + 1 + len(mm.group(1))
        h1 = c2 + 1 + 12
        h2 = h1 + 1 + len(mm.group(3))
        canonical = len(mm.group(1)) == 2 and len(mm.group(3)) == 3
        hits.append(dict(c1=c1, c2=c2, h1=h1, h2=h2, canonical=canonical))
    # drop hits that start inside a previous hit (alternative spacings)
    kept, last_end = [], -1
    for h in sorted(hits, key=lambda x: (x["c1"], not x["canonical"])):
        if h["c1"] > last_end:
            kept.append(h)
            last_end = h["h2"]
    return kept


def runs(zfs):
    """Group ZFs into tandem runs by C1 spacing."""
    if not zfs:
        return []
    out, cur = [], [zfs[0]]
    for a, b in zip(zfs, zfs[1:]):
        if b["c1"] - a["c1"] <= ARRAY_MAX_GAP:
            cur.append(b)
        else:
            out.append(cur)
            cur = [b]
    out.append(cur)
    return out


def identity(a, b):
    n = min(len(a), len(b))
    return sum(x == y for x, y in zip(a[:n], b[:n])) / n if n else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("protein_fasta")
    ap.add_argument("--cds", help="CDS FASTA with identical IDs (codon outputs)")
    ap.add_argument("-o", "--out", default="prdm9")
    ap.add_argument("--min-array", type=int, default=3,
                    help="min tandem ZFs to call the terminal array (default 3)")
    args = ap.parse_args()

    prots = read_fasta(args.protein_fasta)
    cds = read_fasta(args.cds) if args.cds else {}

    f_prot = open(f"{args.out}.per_protein.tsv", "w")
    f_zf = open(f"{args.out}.per_zf.tsv", "w")
    f_ua = open(f"{args.out}.zf_units.faa", "w")
    f_nz = open(f"{args.out}.nonzf.faa", "w")
    f_bc = open(f"{args.out}.binding_codes.tsv", "w")
    f_ar = open(f"{args.out}.array.faa", "w")
    if cds:
        f_nzn = open(f"{args.out}.nonzf.fna", "w")
        f_un = open(f"{args.out}.zf_units.fna", "w")
        f_cc = open(f"{args.out}.zf_contact.fna", "w")
        f_nc = open(f"{args.out}.zf_noncontact.fna", "w")

    f_prot.write("id\tspecies\tprot_len\tn_zf_total\tn_zf_array\tn_noncanonical_array\t"
                 "array_start\tarray_end\tn_zf_upstream\tzf0_present\t"
                 "mean_pairwise_id_array\tbinding_code_string\tcds_ok\n")
    f_zf.write("id\tspecies\tzf_index\tin_array\tarray_pos\tc1\th2\tcanonical\t"
               "unit_aa\tpos-1\tpos+2\tpos+3\tpos+6\tcontact4\n")
    f_bc.write("id\tspecies\tarray_pos\tcontact4\n")

    for pid, prot in prots.items():
        sp = pid.split("|")[0]
        zfs = find_zfs(prot)
        rs = runs(zfs)
        # PRDM9 array = the C-terminal-most run of >= min_array tandem fingers
        array = next((r for r in reversed(rs) if len(r) >= args.min_array), [])
        array_ids = {id(z) for z in array}
        upstream = [z for z in zfs if array and z["h2"] < array[0]["c1"]]

        nt = cds.get(pid, "")
        cds_ok = bool(nt) and len(nt) >= 3 * len(prot) and len(nt) % 3 == 0
        if nt and not cds_ok:
            print(f"WARNING {pid}: CDS length {len(nt)} != 3*protein {len(prot)}; "
                  "codon outputs skipped", file=sys.stderr)

        units, codes = [], []
        for i, z in enumerate(zfs, 1):
            in_arr = id(z) in array_ids
            apos = array.index(z) + 1 if in_arr else ""
            unit = prot[z["c1"]:z["c1"] + UNIT_LEN]
            res = {k: prot[z["h1"] + off] for k, off in CONTACT.items()}
            code = res["-1"] + res["+2"] + res["+3"] + res["+6"]
            f_zf.write(f"{pid}\t{sp}\t{i}\t{int(in_arr)}\t{apos}\t{z['c1']+1}\t"
                       f"{z['h2']+1}\t{int(z['canonical'])}\t{unit}\t"
                       f"{res['-1']}\t{res['+2']}\t{res['+3']}\t{res['+6']}\t{code}\n")
            if not in_arr:
                continue
            units.append(unit)
            codes.append(code)
            tag = f"{pid}|ZF{apos:02d}"
            f_ua.write(f">{tag}\n{unit}\n")
            f_bc.write(f"{pid}\t{sp}\t{apos}\t{code}\n")
            if cds_ok:
                s = z["c1"] * 3
                f_un.write(f">{tag}\n{nt[s:s + 3 * UNIT_LEN]}\n")
                cpos = sorted(z["h1"] + off for off in CONTACT.values())
                upos = range(z["c1"], min(z["c1"] + UNIT_LEN, len(prot)))
                f_cc.write(f">{tag}\n" + "".join(nt[3*p:3*p+3] for p in cpos) + "\n")
                f_nc.write(f">{tag}\n" + "".join(nt[3*p:3*p+3] for p in upos
                                                 if p not in cpos) + "\n")

        pair = [identity(a, b) for a, b in itertools.combinations(units, 2)]
        mean_id = f"{sum(pair)/len(pair):.3f}" if pair else "NA"
        a0, a1 = (array[0]["c1"] + 1, array[-1]["h2"] + 1) if array else ("NA", "NA")
        f_prot.write(f"{pid}\t{sp}\t{len(prot)}\t{len(zfs)}\t{len(array)}\t"
                     f"{sum(not z['canonical'] for z in array)}\t{a0}\t{a1}\t"
                     f"{len(upstream)}\t{int(bool(upstream))}\t{mean_id}\t"
                     f"{'-'.join(codes) or 'NA'}\t{int(cds_ok)}\n")
        nonzf = prot[:array[0]["c1"]] if array else prot
        f_nz.write(f">{pid}\n{nonzf}\n")
        if array:
            f_ar.write(f">{pid}\n{prot[array[0]['c1']:array[-1]['h2'] + 1]}\n")
        if cds_ok:
            f_nzn.write(f">{pid}\n{nt[:3 * len(nonzf)]}\n")

    for f in [f_prot, f_zf, f_ua, f_nz, f_bc, f_ar] + ([f_un, f_cc, f_nc, f_nzn] if cds else []):
        f.close()


if __name__ == "__main__":
    main()
