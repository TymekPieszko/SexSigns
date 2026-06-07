#!/usr/bin/env bash
set -euo pipefail

VCF1=${snakemake_input[vcf1]}
VCF2=${snakemake_input[vcf2]}
VCF_2HAPLOIDS=${snakemake_output[vcf_2haploids]}
VCF_DIPLOID=${snakemake_output[vcf_diploid]}


bcftools merge $VCF1 $VCF2 --missing-to-ref --merge all -Ov | \
bcftools +fixploidy - -- -f 1 > $VCF_2HAPLOIDS && \
awk 'BEGIN{OFS="\t"} /^##/{print;next} /^#CHROM/{sub(/.HP1/,"",$10);NF=10;print;next} {$10=$10"|" $11;NF=10;print}' $VCF_2HAPLOIDS | \
bgzip -c > $VCF_DIPLOID && \
bcftools index $VCF_DIPLOID