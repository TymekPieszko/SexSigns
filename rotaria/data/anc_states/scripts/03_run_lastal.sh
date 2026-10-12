#!/usr/bin/env bash
set -euo pipefail

### Run as:
# ./03_run_lastal.sh R_magnacalcarata_RM9 ../genomes/R_socialis_RS1.fa.gz 18

TARGET=$1
QUERY=$2
THREADS=$3
NAME=$(basename "$QUERY" .fa.gz)
DB_PREFIX=../01_db/${TARGET}/${TARGET}
TRAIN=../02_train/${TARGET}_${NAME}.train
MAF=../03_maf/${TARGET}~${NAME}.maf.gz
mkdir -p $(dirname $MAF)

# Adding --split-f=MAF+ is almost (but not fully) quivalent to a double last-split pass 
# (first time without -r).
# Double split is needed for a one-to-one alignment (single split gives many-to-one).
# Decided not to use maf-linked (but see https://gitlab.com/mcfrith/last/-/blob/main/doc/last-cookbook.rst#aligning-human--chimp-genomes)
lastal -P $THREADS -p $TRAIN -m50 -C2 --split-f=MAF+ $DB_PREFIX $QUERY \
    | last-split -r \
    | gzip -c > $MAF.part
mv $MAF.part $MAF
