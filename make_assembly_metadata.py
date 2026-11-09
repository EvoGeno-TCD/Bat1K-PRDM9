#!/usr/bin/env python3
"""
Build bat_assemblies.tsv from the TOGA2 species table.

Columns: short, assembly, species, common, family, accession, asm_name,
         tech, haplotype, individual, url

tech (heuristic from the accession/assembly name; edit by hand where you know better):
  long_read_phased  assembly name says hap1/hap2/pat/mat
  long_read         ToLID-style name (mMyoMyo1...) or VGP/DToL-like
  short_read        Zoonomia (_BIUU) or old GCF draft (e.g. myoBra1)
  short_read_HiC    DNA Zoo (short-read draft + Hi-C scaffolding)
  unknown           everything else (Bat1K / TBG / Dryad / CNGB etc. - check!)

individual: ToLID stem (mMyoNat1) when present, else the short name with the
            haplotype letter stripped. Assemblies sharing an individual are
            the same animal (haplotypes and/or assembly versions).
"""
import re, sys

TOLID = re.compile(r"(?<![A-Za-z])m[A-Z][a-z]{2}[A-Z][a-z]{2}\d*")
OLD_UCSC = {"myoBra1", "myoDav1", "pteAle1"}
# Bat1K reference genomes (Jebb et al. 2020): PacBio CLR + Hi-C
BAT1K_2020 = {"myoMyo6", "pipKuh2", "molMol2", "phyDis3", "rhiFer5", "rouAeg4"}


def tech_of(short, acc, name, src):
    n = name.lower()
    if "biuu" in n:
        return "short_read"
    if "dna zoo" in src.lower():
        return "short_read_HiC"
    if re.search(r"hap\s*[12]|haplotype\s*[12]|\.pat\b|\.mat\b|_mat_|_pat", n):
        return "long_read_phased"
    if TOLID.search(name) or re.search(r"t2t", n):
        return "long_read"
    if short in BAT1K_2020:
        return "long_read"
    if short in OLD_UCSC:
        return "short_read"
    return "unknown"


def hap_of(short, name):
    n = name.lower()
    m = re.search(r"hap\s*([12])|haplotype\s*([12])", n)
    if m:
        return "hap" + (m.group(1) or m.group(2))
    if re.search(r"\.pat|_pat|paternal", n):
        return "pat"
    if re.search(r"\.mat|_mat|maternal", n):
        return "mat"
    if re.search(r"\d[AB]$", short):
        return "letter" + short[-1]          # haplotype letter, mapping unknown
    return "primary"


def main(table, out, clade="Chiroptera"):
    rows = []
    for line in open(table, encoding="utf-8"):
        f = [x.strip() for x in line.rstrip("\n").split("\t")]
        if len(f) < 4 or clade not in f[3]:
            continue
        sp_common = f[0]
        species = sp_common.split(" (")[0].strip()
        common = sp_common[sp_common.find("(") + 1:sp_common.rfind(")")] if "(" in sp_common else ""
        asm = f[1].split(" (")[0].strip()
        inner = f[1][f[1].find("(") + 1:f[1].rfind(")")]
        acc, _, name = inner.partition(",")
        acc, name = acc.strip(), name.strip()
        short = asm[2:] if asm.startswith("HL") else asm
        fam = next((t.strip() for t in reversed(f[3].split(";")) if t.strip().endswith("idae")), "NA")
        tech = tech_of(short, acc, name, acc)
        hap = hap_of(short, name)
        tol = TOLID.search(name)
        indiv = tol.group(0) if tol else re.sub(r"[AB]$", "", short)
        url = f"http://genome.senckenberg.de/data/{asm}/HLTOGA2integratedv1.bb"
        rows.append([short, asm, species, common, fam, acc or "NA", name or "NA",
                     tech, hap, indiv, url])
    with open(out, "w") as fh:
        fh.write("short\tassembly\tspecies\tcommon\tfamily\taccession\tasm_name\t"
                 "tech\thaplotype\tindividual\turl\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    print(f"{len(rows)} assemblies, {len({r[2] for r in rows})} species -> {out}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("usage: make_assembly_metadata.py <toga_table.txt> <out.tsv> [clade]")
    main(*sys.argv[1:])
