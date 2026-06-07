#!/usr/bin/env bash
set -euo pipefail

# CLADES_FILE=/data/biol-bdelloids/scro4331/SexSigns_2025/yeast/1.args/sticcs/clades_PH.txt

CLADES=(
    12.West_African_cocoa
    3.Brazilian_bioethanol
    '1.Wine_European_(subclade_3)'
    25.Sake
)

# while IFS= read -r CLADE || [[ -n "$CLADE" ]]; do
#     echo $CLADE
#     python calc_ld.py "$CLADE" 1000 20000 100 100
# done < "$CLADES_FILE"

for CLADE in "${CLADES[@]}"; do
    echo "$CLADE"
    python calc_ld.py "$CLADE" 3000 20000 100 100
done