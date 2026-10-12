#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C # Deterministic glob order, identical for R1 and R2

IN_DIR=../0.raw
OUT_DIR=../1.merged

for lib in $IN_DIR/*/*/; do
    lib=${lib%/}
    sp=$(basename $(dirname $lib))
    name=$(basename $lib)
    mkdir -p $OUT_DIR/$sp
    cat $lib/*_1.fastq.gz > $OUT_DIR/$sp/${name}_R1.fastq.gz
    cat $lib/*_2.fastq.gz > $OUT_DIR/$sp/${name}_R2.fastq.gz
    echo $sp/$name
done