#!/usr/bin/env python3

import sys
import csv
import re
from html import unescape
from bs4 import BeautifulSoup

infile = sys.argv[1]
outfile = sys.argv[2]

def clean_html(x):
    if not x:
        return ""
    x = unescape(x)
    soup = BeautifulSoup(x, "html.parser")
    text = soup.get_text(" ", strip=True)
    text = re.sub(r"\s+", " ", text)
    return text

def parse_mutations(html):
    if "<tr>" not in html:
        return ""

    soup = BeautifulSoup(unescape(html), "html.parser")
    muts = []

    for tr in soup.find_all("tr"):
        cells = [td.get_text(strip=True) for td in tr.find_all("td")]
        if len(cells) >= 7:
            exon, pos, mut_type, change, is_inactivating, note, mut_id = cells[:7]
            muts.append(
                f"{mut_id}:{mut_type}:exon{exon}:pos{pos}:{change}:inactivating={is_inactivating}:note={note}"
            )

    return ";".join(muts)

with open(infile) as f, open(outfile, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "query_species_or_db",
        "chrom",
        "start",
        "end",
        "name",
        "strand",
        "query_transcript",
        "gene",
        "projection_id",
        "reference_locus",
        "query_locus",
        "orthology_score",
        "chain_score",
        "nucleotide_divergence",
        "status",
        "perc_intact_ign_M",
        "perc_intact_int_M",
        "intact_codon_prop",
        "mid_intact",
        "mid_pres",
        "inactivating_mutations",
        "protein_alignment_clean"
    ])

    for line in f:
        line = line.rstrip("\n")
        if not line:
            continue

        fields = line.split("\t")

        chrom = fields[0]
        start = fields[1]
        end = fields[2]
        name = fields[3]
        strand = fields[5]

        # name examples:
        # hg38.ENST00000448696.4#IFNE#2
        # HLbosTau10.XM_059889365#IFNE#6
        name_parts = name.split("#")
        gene = name_parts[1] if len(name_parts) > 1 else ""
        projection_id = name_parts[2] if len(name_parts) > 2 else ""

        query_transcript = fields[12] if len(fields) > 12 else ""
        reference_locus = fields[13] if len(fields) > 13 else ""
        query_locus = fields[14] if len(fields) > 14 else ""

        orthology_score = fields[15] if len(fields) > 15 else ""
        chain_score = fields[16] if len(fields) > 16 else ""
        nucleotide_divergence = fields[17] if len(fields) > 17 else ""

        status = fields[22] if len(fields) > 22 else ""

        perc_intact_ign_M = fields[23] if len(fields) > 23 else ""
        perc_intact_int_M = fields[24] if len(fields) > 24 else ""
        intact_codon_prop = fields[25] if len(fields) > 25 else ""
        mid_intact = fields[27] if len(fields) > 27 else ""
        mid_pres = fields[28] if len(fields) > 28 else ""

        protein_alignment_html = fields[29] if len(fields) > 29 else ""
        svg_html = fields[30] if len(fields) > 30 else ""
        transcript_link = fields[31] if len(fields) > 31 else ""
        mutations_html = fields[32] if len(fields) > 32 else ""

        writer.writerow([
            "",  # fill later if looping across species
            chrom,
            start,
            end,
            name,
            strand,
            query_transcript,
            gene,
            projection_id,
            reference_locus,
            query_locus,
            orthology_score,
            chain_score,
            nucleotide_divergence,
            status,
            perc_intact_ign_M,
            perc_intact_int_M,
            intact_codon_prop,
            mid_intact,
            mid_pres,
            parse_mutations(mutations_html),
            clean_html(protein_alignment_html),
        ])
