#!/usr/bin/env bash
# Split reference contigs into N chunks of roughly equal total length.
# Contigs are kept whole and in .fai order, so chunk order == coordinate order
# (required by `bcftools concat --naive`). One chunk per line, contigs
# comma-separated, ready to pass to `bcftools mpileup -r`.
#
# Usage: ROOT=/path/to/calls bash scripts/make_chunks.sh <ref.fa.fai> <n_chunks> > Rmag_chunks.txt
set -euo pipefail
export LC_ALL=C

FAI=$1
N=${2:-50}

TOTAL=$(awk '{s+=$2} END{print s}' $FAI)

awk -v total=$TOTAL -v n=$N '
BEGIN { target = total / n; cur = 0; line = "" }
{
    line = (line == "" ? $1 : line "," $1)
    cur += $2
    if (cur >= target) { print line; line = ""; cur = 0 }
}
END { if (line != "") print line }
' $FAI

