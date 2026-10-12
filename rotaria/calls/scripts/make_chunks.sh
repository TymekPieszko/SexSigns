#!/usr/bin/env bash
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

