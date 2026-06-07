#!/usr/bin/env bash
set -euo pipefail

VCFS_IN=${snakemake_input[vcfs_in]} # Don't double-quote; list needs to be split.
VCF_OUTGROUP=${snakemake_input[vcf_outgroup]}
VCF_OUT=${snakemake_output[vcf_out]}

TMP=${VCF_OUT/.polar/}

bcftools merge $VCFS_IN --missing-to-ref --merge all -Ou | bcftools view - -m2 -M2 -v snps -Oz > $TMP &&
bcftools index $TMP &&
python scripts/polarise.py $TMP $VCF_OUTGROUP $VCF_OUT &&
rm $TMP $TMP.csi