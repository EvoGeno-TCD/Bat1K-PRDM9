#!/usr/bin/env python3
"""
PRDM9 QC and bookkeeping across assemblies.

Inputs
  --meta       bat_assemblies.tsv (make_assembly_metadata.py)
  --rawdir     prdm9_toga_raw/ (run_toga_prdm9.sh: *.PRDM9.bed, *.neighbours.bed, *.none)
  --status-col 1-based column of the TOGA status in the bigBed rows (check bigBedInfo -as)
  --min-chain-score  projections below this chain orthology score (column 16), or with an
               extra flag in the name (e.g. #retro), are counted separately as
               'non_orthologous' and never decide a species' status (default 0.5,
               same as MIN_CHAIN in run_toga_prdm9.sh)
  --zf         <prefix>.per_protein.tsv from prdm9_zf_analysis.py      (optional)
  --domains    prdm9.domains.tsv from prdm9_downstream.sh              (optional)
  --reannot    prdm9_reannot.report.tsv from reannotate_prdm9.py       (optional)
  --cleandir   prdm9_toga_clean/ (parser output with inactivating_mutations) (optional)
  --zf-exon    reference exon number holding the ZF array; default 'auto' =
               the most common exon (block) count among Intact projections,
               i.e. the last coding exon. Printed to stderr - check it.

Outputs (prefix -o)
  .status_by_assembly.tsv  every TOGA status seen per assembly, incl. absent,
                           plus the mutation evidence of its best projection
  .status_by_species.tsv   best status per species, split by assembly tech,
                           with loss evidence and an interpretation
  .mutations.tsv           every projection: inactivating mutations split
                           into ZF exon vs other exons
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
    if s == "non_orthologous":
        return len(STATUS_RANK) + 2
    return STATUS_RANK.index(s) if s in STATUS_RANK else len(STATUS_RANK)


MUT_RE = re.compile(r"^(?P<id>[^:]*):(?P<type>[^:]*):exon(?P<exon>[^:]*):pos(?P<pos>[^:]*):"
                    r"(?P<change>.*?):inactivating=(?P<inact>[^:]*):note=(?P<note>.*)$")
# evidence classes, most to least severe
EVIDENCE = ["outside_zf_exon", "deleted_exons", "missing_exons_only",
            "zf_exon_only", "non_inactivating_only", "none", "no_data"]
EVIDENCE_NOTE = {
    "outside_zf_exon": "inactivating mutation(s) outside the ZF exon: strongest loss evidence",
    "deleted_exons": "exon(s) deleted outside the ZF exon: loss evidence, check alignment chains",
    "missing_exons_only": "only missing exons (often assembly gaps/chain breaks): weak evidence",
    "zf_exon_only": "mutations only in the ZF exon: often a projection artefact - check reannotation",
    "non_inactivating_only": "no inactivating mutations",
    "none": "no mutations reported",
    "no_data": "no parsed projection data",
    "absent": "no PRDM9 projection in TOGA: search the genome directly "
              "(miniprot/tblastn with a bat PRDM9) before calling a loss",
}


def is_true(x):
    return x.strip().lower() in ("true", "yes", "1", "+", "t", "y")


def parse_mutations(field, zf_exon):
    """-> dict summary of one projection's inactivating_mutations string."""
    muts = []
    for m in filter(None, (x.strip() for x in (field or "").split(";"))):
        mm = MUT_RE.match(m)
        if mm:
            muts.append(mm.groupdict())
    inact = [m for m in muts if is_true(m["inact"])]
    def in_zf(m):
        return zf_exon is not None and m["exon"].strip() == str(zf_exon)
    zf = [m for m in inact if in_zf(m)]
    out = [m for m in inact if not in_zf(m)]
    def kind(m):
        t = (m["type"] + " " + m["change"]).lower()
        if "missing" in t and "exon" in t:
            return "missing"
        if "deleted" in t or ("del" in t and "exon" in t and "fs" not in t):
            return "deleted"
        return "other"
    if not muts:
        ev = "none"
    elif not inact:
        ev = "non_inactivating_only"
    elif any(kind(m) == "other" for m in out):
        ev = "outside_zf_exon"
    elif any(kind(m) == "deleted" for m in out):
        ev = "deleted_exons"
    elif out:
        ev = "missing_exons_only"
    else:
        ev = "zf_exon_only"
    return dict(n_mut=len(muts), n_inact=len(inact), n_zf=len(zf), n_out=len(out),
                exons_out=",".join(sorted({m["exon"] for m in out}, key=lambda x: (len(x), x))),
                types_out=",".join(sorted({m["type"] for m in out})),
                types_zf=",".join(sorted({m["type"] for m in zf})), evidence=ev)


def is_orth(name, chain, min_chain):
    """orthologous projection: no extra flag after tx#GENE#id, chain score >= min."""
    if len([p for p in name.split("#")]) > 3:
        return False
    try:
        return float(chain) >= min_chain
    except ValueError:
        return True          # no score available: keep


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
    ap.add_argument("--cleandir")
    ap.add_argument("--min-chain-score", type=float, default=0.5)
    ap.add_argument("--zf-exon", default="auto")
    ap.add_argument("--n-flank-genes", type=int, default=5)
    ap.add_argument("-o", "--out", default="prdm9_qc")
    a = ap.parse_args()
    sc = a.status_col - 1

    meta = {r["short"]: r for r in read_tsv(a.meta)}

    # ---------------- 1. status, including losses / absence -----------------
    by_asm, blocks_intact = {}, Counter()
    for short in meta:
        bed = os.path.join(a.rawdir, f"{short}.PRDM9.bed")
        if os.path.exists(os.path.join(a.rawdir, f"{short}.none")):
            by_asm[short] = Counter({"absent_from_TOGA": 1})
        elif os.path.exists(bed):
            c = Counter()
            for line in open(bed):
                f = line.rstrip("\n").split("\t")
                if not is_orth(f[3] if len(f) > 3 else "", f[15] if len(f) > 15 else "",
                               a.min_chain_score):
                    c["non_orthologous"] += 1
                    continue
                st = f[sc].strip() if len(f) > sc else "NA"
                c[st] += 1
                if st.lower() == "intact" and len(f) > 9 and f[9].isdigit():
                    blocks_intact[int(f[9])] += 1
            by_asm[short] = c
        else:
            by_asm[short] = Counter({"not_processed": 1})

    # ---------------- 1b. inactivating mutations (parser output) -------------
    if a.zf_exon == "auto":
        zf_exon = blocks_intact.most_common(1)[0][0] if blocks_intact else None
        print(f"ZF exon (auto, most common exon count of Intact projections): {zf_exon} "
              f"[counts: {dict(blocks_intact.most_common(5))}] - override with --zf-exon",
              file=sys.stderr)
    else:
        zf_exon = int(a.zf_exon)
    if zf_exon is None and a.cleandir:
        print("WARNING: ZF exon unknown (no Intact projections with an exon count); every "
              "inactivating mutation will count as 'outside ZF exon'. Set --zf-exon.",
              file=sys.stderr)
    csv.field_size_limit(sys.maxsize)
    best_proj = {}            # short -> mutation summary of best projection
    if a.cleandir:
        with open(a.out + ".mutations.tsv", "w") as mf:
            mf.write("short\tspecies\tname\tstatus\tn_mutations\tn_inactivating\t"
                     "n_inact_zf_exon\tn_inact_other_exons\tother_exons_hit\t"
                     "types_other_exons\ttypes_zf_exon\tevidence\n")
            for short in meta:
                path = os.path.join(a.cleandir, f"{short}.PRDM9.clean.tsv")
                if not os.path.exists(path):
                    continue
                cands = []
                with open(path) as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        if not is_orth(r.get("name", ""), r.get("chain_score", "1"),
                                       a.min_chain_score):
                            continue
                        ms = parse_mutations(r.get("inactivating_mutations", ""), zf_exon)
                        st = r.get("status", "NA")
                        cands.append((status_rank(st), EVIDENCE.index(ms["evidence"]) * -1,
                                      ms["n_inact"], st, r.get("name", "NA"), ms))
                        mf.write(f"{short}\t{meta[short]['species']}\t{r.get('name','NA')}\t{st}\t"
                                 f"{ms['n_mut']}\t{ms['n_inact']}\t{ms['n_zf']}\t{ms['n_out']}\t"
                                 f"{ms['exons_out'] or '-'}\t{ms['types_out'] or '-'}\t"
                                 f"{ms['types_zf'] or '-'}\t{ms['evidence']}\n")
                if cands:
                    # best = best status, then the least severe evidence, then fewest mutations
                    best = min(cands, key=lambda x: (x[0], x[1], x[2]))
                    best_proj[short] = dict(best[5], status=best[3], name=best[4])

    def evidence_of(short):
        if short in best_proj:
            return best_proj[short]["evidence"]
        return "no_data" if not by_asm[short].get("absent_from_TOGA") else "absent"

    with open(a.out + ".status_by_assembly.tsv", "w") as fh:
        fh.write("short\tspecies\tfamily\ttech\tbest_status\tall_statuses\t"
                 "best_projection_evidence\tother_exons_hit\ttypes_other_exons\ttypes_zf_exon\n")
        for s_, c in by_asm.items():
            best = min(c, key=status_rank)
            m = meta[s_]
            bp = best_proj.get(s_, {})
            fh.write(f"{s_}\t{m['species']}\t{m['family']}\t{m['tech']}\t{best}\t"
                     + ";".join(f"{k}:{v}" for k, v in c.most_common())
                     + f"\t{evidence_of(s_)}\t{bp.get('exons_out') or '-'}\t"
                     f"{bp.get('types_out') or '-'}\t{bp.get('types_zf') or '-'}\n")

    sp_best = defaultdict(dict)       # species -> grp -> (status, short)
    for s_, c in by_asm.items():
        m = meta[s_]
        grp = "long" if m["tech"].startswith("long") else "other"
        best = min(c, key=status_rank)
        prev = sp_best[m["species"]].get(grp)
        ev = evidence_of(s_)
        evr = EVIDENCE.index(ev) if ev in EVIDENCE else len(EVIDENCE)
        # prefer better status; at equal status prefer the LEAST severe evidence
        key = (status_rank(best), -evr)
        if prev is None or key < prev[0]:
            sp_best[m["species"]][grp] = (key, best, s_, ev)
    with open(a.out + ".status_by_species.tsv", "w") as fh:
        fh.write("species\tfamily\tbest_status_longread\tbest_status_other\tbest_overall\t"
                 "evidence_assembly\tloss_evidence\tinterpretation\tnote\n")
        fam = {m["species"]: m["family"] for m in meta.values()}
        for sp, d in sorted(sp_best.items(), key=lambda x: (fam[x[0]], x[0])):
            pick = d.get("long") or d.get("other")        # evidence from long-read if any
            allb = min((v[1] for v in d.values()), key=status_rank)
            note = []
            if "long" in d and "other" in d and status_rank(d["other"][1]) < status_rank(d["long"][1]):
                note.append("short-read better than long-read: check")
            if status_rank(allb) > 1:
                note.append("candidate loss/truncation")
            ev = pick[3]
            fh.write(f"{sp}\t{fam[sp]}\t{d['long'][1] if 'long' in d else 'NA'}\t"
                     f"{d['other'][1] if 'other' in d else 'NA'}\t{allb}\t{pick[2]}\t{ev}\t"
                     f"{EVIDENCE_NOTE.get(ev, ev)}\t{'; '.join(note)}\n")

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
