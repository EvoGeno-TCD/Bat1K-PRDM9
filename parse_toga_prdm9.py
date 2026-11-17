#!/usr/bin/env python3
"""
Turn TOGA2 bigBed rows (38-column 'togaBigBed' schema) into a readable table.
Replaces the lab's parse_toga_ifne.py: standard library only (no bs4), and the
column names follow the bigBed schema (bigBedInfo -as).

  parse_toga_prdm9.py rows.bed out.tsv

Output columns (first column is filled in later by run_toga_prdm9.sh):
  query_species_or_db chrom start end name strand ref_trans_id gene projection_id
  flag ref_region query_region chain_score chain_synteny chain_flank
  chain_gl_cds_fract chain_loc_cds_fract chain_exon_cov chain_intron_cov status
  perc_intact_ign_M perc_intact_int_M intact_codon_prop ouf_prop mid_intact mid_pres
  num_affected_exon perc_affected_exons inactivating_mutations cds protein

inactivating_mutations: ';'-separated records
  <mut_id>:<type>:exon<exon>:pos<pos>:<change>:inactivating=<YES/NO>:note=<note>
"""
import html, re, sys

TR = re.compile(r"<tr>(.*?)</tr>", re.S | re.I)
TD = re.compile(r"<td>(.*?)</td>", re.S | re.I)
TAG = re.compile(r"<[^>]+>")

OUT = ["query_species_or_db", "chrom", "start", "end", "name", "strand", "ref_trans_id",
       "gene", "projection_id", "flag", "ref_region", "query_region", "chain_score",
       "chain_synteny", "chain_flank", "chain_gl_cds_fract", "chain_loc_cds_fract",
       "chain_exon_cov", "chain_intron_cov", "status", "perc_intact_ign_M",
       "perc_intact_int_M", "intact_codon_prop", "ouf_prop", "mid_intact", "mid_pres",
       "num_affected_exon", "perc_affected_exons", "inactivating_mutations", "cds", "protein"]


def cell(x):
    return html.unescape(TAG.sub("", x)).replace("\xa0", " ").replace("\t", " ").strip()


def mutations(htm):
    out = []
    for tr in TR.findall(htm or ""):
        c = [cell(x) for x in TD.findall(tr)]
        if len(c) >= 7:
            exon, pos, typ, change, inact, note, mid = c[:7]
            out.append(f"{mid}:{typ}:exon{exon}:pos{pos}:{change}:inactivating={inact}:note={note}")
    return ";".join(out)


def main(inp, outp):
    n = bad = 0
    with open(inp) as fh, open(outp, "w") as out:
        out.write("\t".join(OUT) + "\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 38:
                bad += 1
                continue
            parts = f[3].split("#")
            gene = parts[1] if len(parts) > 1 else ""
            proj = parts[2] if len(parts) > 2 else ""
            flag = "#".join(parts[3:]) if len(parts) > 3 else ""   # e.g. retro
            row = ["", f[0], f[1], f[2], f[3], f[5], f[12], gene, proj, flag, f[13], f[14],
                   f[15], f[16], f[17], f[18], f[19], f[20], f[21], f[22], f[23], f[24],
                   f[25], f[26], f[27], f[28], f[36], f[37], mutations(f[32]),
                   f[34].strip(), f[35].strip()]
            out.write("\t".join(row) + "\n")
            n += 1
    print(f"parsed {n} rows" + (f", skipped {bad} with <38 columns" if bad else ""),
          file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
