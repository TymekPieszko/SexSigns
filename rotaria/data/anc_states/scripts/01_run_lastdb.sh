#!/usr/bin/env bash
set -euo pipefail

### Run as:
# ./01_run_lastdb.sh R_magnacalcarata_RM9

TARGET=$1
REF=../genomes/${TARGET}.fa.gz
DB_PREFIX=../01_db/${TARGET}/${TARGET}
THREADS=9

mkdir -p $(dirname $DB_PREFIX)

# https://gitlab.com/mcfrith/last/-/blob/main/doc/lastdb.rst
# -uMAM8 : Seeding scheme that finds weak similarities with high sensitivity (https://gitlab.com/mcfrith/last/-/blob/main/doc/last-seeds.rst)
# -R01 : Find simple repeats with tantan and soft-mask (lowercase) them, keping them in the index so alignments can extend through. THIS IS DEFAULT AND NOT NEEDED.
# -c : Soft-mask lowercase regions during seeding.
# -v : Verbose
lastdb -v -P $THREADS -uMAM8 -R01 -c $DB_PREFIX $REF