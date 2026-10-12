#!/usr/bin/env bash
set -euo pipefail

### Run as:
# TARGET=R_magnacalcarata_RM9
# QUERY=../genomes/R_socialis_RS1
# THREADS=9
# ./02_run_last-train.sh $TARGET $QUERY $THREADS

TARGET=$1
QUERY=$2
THREADS=$3
DB_PREFIX=../01_db/${TARGET}/${TARGET}
OUT=../02_train/${TARGET}_$(basename $QUERY .fa.gz).train
mkdir -p $(dirname $OUT)

# https://gitlab.com/mcfrith/last/-/blob/main/doc/last-train.rst
# Flags after https://github.com/mcfrith/last-genome-alignments#2017-human-ape-alignments
# The three 'sym' flags enforce substitution / indel symmetry.
last-train -P $THREADS --revsym --matsym --gapsym -E 0.05 -C 2 $DB_PREFIX $QUERY > $OUT
# Alternative for 'more distantly related genomes':
# https://github.com/mcfrith/last