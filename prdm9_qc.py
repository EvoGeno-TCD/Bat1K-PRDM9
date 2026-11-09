#!/usr/bin/env python3
"""
PRDM9 QC and bookkeeping across assemblies.

Inputs
  --meta       bat_assemblies.tsv (make_assembly_metadata.py)
  --rawdir     prdm9_toga_raw/ (run_toga_prdm9.sh: *.PRDM9.bed, *.neighbours.bed, *.none)
  --status-col 1-based column of the TOGA status in the bigBed rows (check bigBedInfo -as)
  --zf         <prefix>.per_protein.tsv from prdm9_zf_analysis.py      (optional)
  --domains    prdm9.domains.tsv from prdm9_downstream.sh              (optional)
  --reannot    prdm9_reannot.report.tsv from reannotate_prdm9.py       (optional)

Outputs (prefix -o)
  .status_by_assembly.tsv  every TOGA status seen per assembly, incl. absent
  .status_by_species.tsv   best status per species, split by assembly tech
  .synteny.tsv             flanking genes per locus and a synteny score
  .master.tsv              one row per protein: ZF/domains/reannot flags/meta/gate
  .haplotypes.tsv          per individual: hap1 vs hap2 (allelic) and version
                           agreement (assembly QC)
  .representatives.tsv     one protein per species for species-level analyses
"""
import argparse, csv, glob, os, re, sys
from collections import Counter, defaultdict

TECH_RANK = {"long_read_phased": 0, "long_read": 1, "unknown": 2,
             "short_read_HiC": 3, "short_read": 4}
STATUS_RANK = ["intact", "partially intact", "uncertain loss", "missing", "lost",
               "paralogous projection"]


def gene_of(name):
    toks = [t for t in name.split("#") if t]
    return toks[1] if len(toks) > 1 else toks[0]


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def status_rank(s):
    s = s.lower().strip()
    return STATUS_RANK.index(s) if s in STATUS_RANK else len(STATUS_RANK)


def version_of(short):
    m = re.search(r"(\d+)[AB]?$", short)
    return int(m.group(1)) if m else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--rawdir", required=True)
    ap.add_argument("--status-col", type=int, default=23)
    ap.add_argument("--zf")
    ap.add_argument("--domains")
    ap.add_argument("--reannot")
    ap.add_argument("--n-flank-genes", type=int, default=5)
    ap.add_argument("-o", "--out", default="prdm9_qc")
    a = ap.parse_args()
    sc = a.status_col - 1

    meta = {r["short"]: r for r in read_tsv(a.meta)}

    # ---------------- 1. status, including losses / absence -----------------
    by_asm = {}
    for short in meta:
        bed = os.path.join(a.rawdir, f"{short}.PRDM9.bed")
        if os.path.exists(os.path.join(a.rawdir, f"{short}.none")):
            by_asm[short] = Counter({"absent_from_TOGA": 1})
        elif os.path.exists(bed):
            c = Counter()
            for line in open(bed):
                f = line.rstrip("\n").split("\t")
                c[f[sc].strip() if len(f) > sc else "NA"] += 1
            by_asm[short] = c
        else:
            by_asm[short] = Counter({"not_processed": 1})
    with open(a.out + ".status_by_assembly.tsv", "w") as fh:
        fh.write("short\tspecies\tfamily\ttech\tbest_status\tall_statuses\n")
        for s, c in by_asm.items():
            best = min(c, key=status_rank)
            m = meta[s]
            fh.write(f"{s}\t{m['species']}\t{m['family']}\t{m['tech']}\t{best}\t"
                     + ";".join(f"{k}:{v}" for k, v in c.most_common()) + "\n")

    sp_best = defaultdict(dict)
    for s, c in by_asm.items():
        m = meta[s]
        grp = "long" if m["tech"].startswith("long") else "other"
        best = min(c, key=status_rank)
        prev = sp_best[m["species"]].get(grp)
        if prev is None or status_rank(best) < status_rank(prev):
            sp_best[m["species"]][grp] = best
    with open(a.out + ".status_by_species.tsv", "w") as fh:
        fh.write("species\tfamily\tbest_status_longread\tbest_status_other\tbest_overall\tnote\n")
        fam = {m["species"]: m["family"] for m in meta.values()}
        for sp, d in sorted(sp_best.items(), key=lambda x: (fam[x[0]], x[0])):
            allb = min(d.values(), key=status_rank)
            note = ""
            if "long" in d and "other" in d and status_rank(d["other"]) < status_rank(d["long"]):
                note = "short-read better than long-read: check"
            if status_rank(allb) > 1:
                note = (note + "; " if note else "") + "candidate loss/truncation"
            fh.write(f"{sp}\t{fam[sp]}\t{d.get('long','NA')}\t{d.get('other','NA')}\t{allb}\t{note}\n")

    # ---------------- 2. synteny ---------------------------------------------
    loci = {}   # (short, locus) -> dict(left=[...], right=[...])
    for nb in glob.glob(os.path.join(a.rawdir, "*.PRDM9.neighbours.bed")):
        short = os.path.basename(nb).replace(".PRDM9.neighbours.bed", "")
        rows = defaultdict(list)
        for line in open(nb):
            f = line.rstrip("\n").split("\t")
            rows[f[0]].append((int(f[2]), int(f[3]), gene_of(f[4])))
        for locus, rr in rows.items():
            chrom, se = locus.rsplit(":", 1)
            ls, le = map(int, se.split("-"))
            genes = {}
            for s, e, g in rr:
                if "prdm9" in g.lower():
                    continue
                d = (ls - e) if e <= ls else (s - le) if s >= le else 0
                side = "L" if e <= ls else "R" if s >= le else "O"
                if g not in genes or abs(d) < abs(genes[g][0]):
                    genes[g] = (d, side)
            left = [g for g, (d, sd) in sorted(genes.items(), key=lambda x: x[1][0]) if sd == "L"]
            right = [g for g, (d, sd) in sorted(genes.items(), key=lambda x: x[1][0]) if sd == "R"]
            loci[(short, locus)] = dict(L=left[:a.n_flank_genes], R=right[:a.n_flank_genes])
    freq = Counter(g for v in loci.values() for g in set(v["L"] + v["R"]))
    core = {g for g, n in freq.items() if loci and n >= 0.5 * len(loci)}
    syn = {}
    with open(a.out + ".synteny.tsv", "w") as fh:
        fh.write("short\tlocus\tleft_genes\tright_genes\tcore_shared\tcore_size\tsynteny_score\tflag\n")
        for (short, locus), v in sorted(loci.items()):
            shared = len(core & set(v["L"] + v["R"]))
            score = shared / len(core) if core else float("nan")
            flag = "SYNTENIC" if score >= 0.5 else "WEAK" if score >= 0.2 else "NON_SYNTENIC"
            syn[(short, locus)] = (score, flag)
            fh.write(f"{short}\t{locus}\t{','.join(v['L'])}\t{','.join(v['R'])}\t"
                     f"{shared}\t{len(core)}\t{score:.2f}\t{flag}\n")
    print(f"core flanking genes (>=50% of loci): {','.join(sorted(core)) or 'none'}", file=sys.stderr)

    if not a.zf:
        return

    # ---------------- 3. master table + PRDM9 gate ----------------------------
    zf = read_tsv(a.zf)
    dom = {r["id"]: r for r in read_tsv(a.domains)} if a.domains else {}
    rea = {r["id"]: r for r in read_tsv(a.reannot)} if a.reannot else {}

    def locus_syn(short, pid):
        m = re.search(r"\|([^|:]+):(\d+)-(\d+)", pid)
        if not m:
            return ("NA", "NA")
        c, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        for (sh, loc), v in syn.items():
            lc, se = loc.rsplit(":", 1)
            ls, le = map(int, se.split("-"))
            if sh == short and lc == c and s < le + 200000 and e > ls - 200000:
                return v
        return ("NA", "NA")

    master = []
    for r in zf:
        pid = r["id"]
        short = pid.split("|")[0]
        m = meta.get(short, {})
        d = dom.get(pid, {})
        rr = rea.get(pid, {})
        score, sflag = locus_syn(short, pid)
        has_set, has_ssx = d.get("SET", "NA"), d.get("SSXRD", "NA")
        reasons = []
        if d and has_set != "1":
            reasons.append("no_SET")
        if d and has_ssx != "1":
            reasons.append("no_SSXRD")
        if int(r["n_zf_array"]) < 3:
            reasons.append("no_ZF_array")
        if sflag == "NON_SYNTENIC":
            reasons.append("non_syntenic")
        rflag = rr.get("flag", "NA")
        if any(x in rflag for x in ("GAP_IN_ZF_EXON", "NO_STOP", "INTERNAL_STOP")):
            reasons.append("zf_exon_unreliable")
        gate = "PRDM9" if not reasons else ("PRDM9_like" if d.get("SET") == "1" else "not_PRDM9")
        master.append(dict(
            id=pid, short=short, species=m.get("species", "NA"), family=m.get("family", "NA"),
            tech=m.get("tech", "NA"), haplotype=m.get("haplotype", "NA"),
            individual=m.get("individual", "NA"), n_zf_array=r["n_zf_array"],
            n_zf_total=r["n_zf_total"], zf0=r["zf0_present"],
            mean_id_array=r["mean_pairwise_id_array"], binding_codes=r["binding_code_string"],
            KRAB=d.get("KRAB", "NA"), SSXRD=has_ssx, SET=has_set,
            synteny=score if score == "NA" else f"{score:.2f}", synteny_flag=sflag,
            reannot_flag=rflag, gate=gate, gate_reasons=",".join(reasons) or "pass",
            zf_count_trustworthy=int(m.get("tech", "").startswith("long")
                                     and "zf_exon_unreliable" not in reasons)))
    cols = list(master[0].keys()) if master else []
    with open(a.out + ".master.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in master:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # main copy per assembly: gate pass first, then longest array
    main_copy = {}
    for row in master:
        key = row["short"]
        rank = (row["gate"] != "PRDM9", -int(row["n_zf_array"]))
        if key not in main_copy or rank < main_copy[key][0]:
            main_copy[key] = (rank, row)
    main_rows = [v[1] for v in main_copy.values()]
    copies = Counter(r["short"] for r in master if r["gate"] == "PRDM9")

    # ---------------- 4. haplotypes and versions -----------------------------
    by_ind = defaultdict(list)
    for row in main_rows:
        by_ind[(row["species"], row["individual"])].append(row)
    with open(a.out + ".haplotypes.tsv", "w") as fh:
        fh.write("species\tindividual\tassemblies\tzf_counts\tversion_check\t"
                 "allelic_zf_diff\tallelic_codes_identical\n")
        for (sp, ind), rows in sorted(by_ind.items()):
            by_hap = defaultdict(list)
            for r in rows:
                by_hap[r["haplotype"]].append(r)
            vcheck = []
            for h, rr in by_hap.items():
                if len(rr) > 1:
                    same = len({(x["n_zf_array"], x["binding_codes"]) for x in rr}) == 1
                    vcheck.append(f"{h}:{'agree' if same else 'DISAGREE'}")
            haps = {h: sorted(rr, key=lambda x: -version_of(x["short"]))[0]
                    for h, rr in by_hap.items()}
            pair = [haps[h] for h in ("hap1", "hap2", "pat", "mat", "letterA", "letterB")
                    if h in haps]
            if len(pair) >= 2:
                adiff = abs(int(pair[0]["n_zf_array"]) - int(pair[1]["n_zf_array"]))
                acode = int(pair[0]["binding_codes"] == pair[1]["binding_codes"])
            else:
                adiff, acode = "NA", "NA"
            fh.write(f"{sp}\t{ind}\t{','.join(r['short'] for r in rows)}\t"
                     f"{','.join(r['n_zf_array'] for r in rows)}\t"
                     f"{';'.join(vcheck) or 'single_version'}\t{adiff}\t{acode}\n")

    # ---------------- 5. one representative per species -----------------------
    by_sp = defaultdict(list)
    for row in main_rows:
        by_sp[row["species"]].append(row)
    with open(a.out + ".representatives.tsv", "w") as fh:
        fh.write("species\tfamily\tid\tshort\ttech\tn_zf_array\tgate\tn_prdm9_copies\t"
                 "n_assemblies\talternatives\n")
        for sp, rows in sorted(by_sp.items()):
            rows.sort(key=lambda r: (r["gate"] != "PRDM9", r["reannot_flag"] not in ("OK", "NA"),
                                     TECH_RANK.get(r["tech"], 9),
                                     r["haplotype"] not in ("hap1", "primary", "pat"),
                                     -version_of(r["short"])))
            b = rows[0]
            fh.write(f"{sp}\t{b['family']}\t{b['id']}\t{b['short']}\t{b['tech']}\t"
                     f"{b['n_zf_array']}\t{b['gate']}\t{copies.get(b['short'], 0)}\t{len(rows)}\t"
                     f"{','.join(r['short'] for r in rows[1:])}\n")


if __name__ == "__main__":
    main()
