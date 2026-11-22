#!/usr/bin/env bash
# Check that this machine can run the pipeline, before starting a long run.
# Nothing is downloaded except a few small test regions.
#   bash check_setup.sh
ok=0; bad=0
say()  { printf '%-34s %s\n' "$1" "$2"; }
pass() { say "$1" "OK  $2"; ok=$((ok+1)); }
fail() { say "$1" "FAIL $2"; bad=$((bad+1)); }

echo "== programs =="
for t in bigBedToBed bigBedInfo twoBitToFa twoBitInfo curl python3 miniprot mafft iqtree2 hyphy Rscript interproscan.sh; do
    if command -v "$t" >/dev/null 2>&1; then pass "$t" "$(command -v "$t")"
    else case "$t" in
            bigBedToBed|bigBedInfo|twoBitToFa|python3) fail "$t" "REQUIRED - install it" ;;
            *) say "$t" "-   missing (that step will be skipped)" ;;
         esac
    fi
done

echo
echo "== network =="
TB="https://genome.senckenberg.de/download/TOGA2/Genomes/Mammalia/Rhinolophus_ferrumequinum__greater_horseshoe_bat__HLrhiFer5__GCA_014108255.1/HLrhiFer5.2bit"
BB="http://genome.senckenberg.de/data/HLrhiFer5/HLTOGA2integratedv1.bb"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT

if bigBedInfo -udcDir="$tmp/u1" "$BB" >"$tmp/bb.txt" 2>"$tmp/bb.err"; then
    pass "TOGA2 annotations (bigBed)" "$(grep -m1 itemCount "$tmp/bb.txt" | tr -d ' ')"
else
    fail "TOGA2 annotations (bigBed)" "$(head -1 "$tmp/bb.err")"
fi

if twoBitInfo -udcDir="$tmp/u2" "$TB" stdout 2>"$tmp/tb.err" | head -3 >"$tmp/tb.txt" && [[ -s "$tmp/tb.txt" ]]; then
    pass "genome sequence (2bit)" "first scaffold: $(head -1 "$tmp/tb.txt" | tr '\t' ' ')"
else
    fail "genome sequence (2bit)" "$(head -1 "$tmp/tb.err")"
fi

if command -v curl >/dev/null 2>&1; then
    code=$(curl -s -o "$tmp/ncbi.fa" -w '%{http_code}' --max-time 60 \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000001.11&rettype=fasta&seq_start=1&seq_stop=200" || echo 000)
    if [[ "$code" == 200 ]] && head -c1 "$tmp/ncbi.fa" | grep -q '>'; then
        pass "NCBI (fallback only)" "reachable"
    else
        say "NCBI (fallback only)" "-   not reachable (HTTP $code) - fine if the 2bit test passed"
    fi
fi

echo
echo "== files in $(pwd) =="
for f in bat_toga_urls.tsv bat_assemblies.tsv parse_toga_prdm9.py; do
    [[ -s "$f" ]] && pass "$f" "$(wc -l < "$f") lines" || fail "$f" "missing"
done
for f in genome_2bit_urls.tsv query_prdm9.faa; do
    if [[ -s "$f" ]]; then pass "$f" "$(wc -l < "$f") lines"
    else say "$f" "-   missing ($([[ $f == genome* ]] && echo 'run: python3 make_genome_urls.py' || echo 'needed for step 2'))"
    fi
done

echo
echo "$ok checks passed, $bad required checks failed."
[[ $bad -eq 0 ]] || echo "Fix the FAIL lines before running the pipeline."
