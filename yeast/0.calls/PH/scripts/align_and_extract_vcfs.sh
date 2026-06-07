#!/usr/bin/env bash
set -euo pipefail

# Loop version - but you'd need to subtract one from the output index!
# for i in 1 2; do
#     echo "${snakemake_input[0]}"
#     echo "${snakemake_input[$i]}"
#     echo "${snakemake_output[$i]}"
#     minimap2 -cx asm5 -t ${snakemake[threads]} --cs "${snakemake_input[0]}" "${snakemake_input[$i]}" > "${snakemake_output[$i]}"
# done

## Assign variables
THREADS=${snakemake[threads]}
REF=${snakemake_input[ref]}
###
HAP1=${snakemake_input[hap1]}
PAF1=${snakemake_output[paf1]}
VCF1=${snakemake_output[vcf1]}
###
HAP2=${snakemake_input[hap2]}
PAF2=${snakemake_output[paf2]}
VCF2=${snakemake_output[vcf2]}
###

echo $HAP1
echo $PAF1
echo $VCF1
echo $HAP1
echo $HAP1
echo $HAP1

### Haplome 1
minimap2 -cx asm5 -t $THREADS --cs $REF $HAP1 > $PAF1 && \
sort -k6,6 -k8,8n $PAF1 | paftools.js call -f $REF -s ${snakemake_wildcards[SAMPLE]}.HP1 -L 5000 - | bgzip -c > $VCF1 && \
bcftools index $VCF1

### Haplome 2
minimap2 -cx asm5 -t $THREADS --cs $REF $HAP2 > $PAF2 && \
sort -k6,6 -k8,8n $PAF2 | paftools.js call -f $REF -s ${snakemake_wildcards[SAMPLE]}.HP2 -L 5000 - | bgzip -c > $VCF2 && \
bcftools index $VCF2