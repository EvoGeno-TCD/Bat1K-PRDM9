#!/usr/bin/env python3
"""
Build genome_2bit_urls.tsv: <short assembly name> <TAB> <URL of its .2bit>

The TOGA2 annotation files can be found under  /data/<assembly>/HLTOGA2integratedv1.bb
but the genome SEQUENCES are somewhere else, e.g.

  https://genome.senckenberg.de/download/TOGA2/Genomes/Mammalia/
      Rhinolophus_ferrumequinum__greater_horseshoe_bat__HLrhiFer5__GCA_014108255.1/HLrhiFer5.2bit

So we download the directory index once and match each assembly ID to its folder.

  python3 make_genome_urls.py                        # downloads the index
  python3 make_genome_urls.py --index-file idx.html  # use a saved copy
  python3 make_genome_urls.py --verify               # also HEAD each .2bit URL

Reads assembly IDs from bat_toga_urls.tsv (the HL name in each URL).
"""
import argparse, os, re, sys, urllib.request

INDEX = "https://genome.senckenberg.de/download/TOGA2/Genomes/Mammalia/"
HREF = re.compile(r'href="([^"?/][^"]*)/"')


def get(url, timeout=120, head=False):
    req = urllib.request.Request(url, method="HEAD" if head else "GET",
                                 headers={"User-Agent": "prdm9-pipeline/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.status if head else r.read().decode("utf-8", "replace")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--urls", default="bat_toga_urls.tsv")
    ap.add_argument("--index-url", default=INDEX)
    ap.add_argument("--index-file")
    ap.add_argument("--verify", action="store_true",
                    help="check each .2bit URL exists (one HEAD request each)")
    ap.add_argument("-o", "--out", default="genome_2bit_urls.tsv")
    a = ap.parse_args()

    if a.index_file:
        html = open(a.index_file, encoding="utf-8", errors="replace").read()
    else:
        print(f"downloading {a.index_url}", file=sys.stderr)
        try:
            html = get(a.index_url)
        except Exception as e:
            sys.exit(f"could not download the index ({e}).\n"
                     f"Save it by hand and re-run with --index-file:\n"
                     f"  curl -sSf '{a.index_url}' -o idx.html")
    folders = [f for f in HREF.findall(html) if "__" in f]
    print(f"{len(folders)} genome folders in the index", file=sys.stderr)

    # assembly id -> folder.  Folder names end with  __<asm>  or  __<asm>__<accession>
    by_asm = {}
    for f in folders:
        for i, tok in enumerate(f.split("__")):
            if i > 0 and tok:
                by_asm.setdefault(tok, f)

    asms = []
    for line in open(a.urls):
        p = line.split()
        if len(p) < 2 or p[0] in ("short", "#") or line.startswith("#"):
            continue
        asms.append((p[0], os.path.basename(os.path.dirname(p[1]))))

    base = a.index_url if a.index_url.endswith("/") else a.index_url + "/"
    found = missing = bad = 0
    with open(a.out, "w") as out:
        for short, asm in asms:
            folder = by_asm.get(asm)
            if folder is None:
                print(f"NO FOLDER for {asm} ({short})", file=sys.stderr)
                missing += 1
                continue
            url = f"{base}{folder}/{asm}.2bit"
            if a.verify:
                try:
                    get(url, timeout=60, head=True)
                except Exception as e:
                    print(f"MISSING FILE {asm}: {url} ({e})", file=sys.stderr)
                    bad += 1
                    continue
            out.write(f"{short}\t{url}\n")
            found += 1
    print(f"wrote {a.out}: {found} assemblies"
          + (f", {missing} with no folder" if missing else "")
          + (f", {bad} whose .2bit was not found" if bad else ""), file=sys.stderr)
    if missing or bad:
        print("Assemblies not listed here fall back to NCBI in run_toga_prdm9.sh.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
