#!/usr/bin/env bash
set -euo pipefail

### Run as:
# ./04_polarise.sh Rmag Rsoc
# ./04_polarise.sh Rsoc Rmag

### Possible simplifications / efficiency fixes FOR THE FUTURE:
# Map only the SNP positions, not every aligned base (the map_sites.py step).
# Can use `maf-convert chain` and `CrossMap vcf`.
# `bcftools isec`/`annotate` could be possibly useful in make_aa.py.

TARGET=$1
COMPARISON=$2
VCF1=../../vcfs_total/${TARGET}_phased.vcf.gz
VCF2=../../vcfs_total/${COMPARISON}_phased.vcf.gz
MAF=../03_maf/R_magnacalcarata_RM9~R_socialis_RS1.maf.gz
OUT_DIR=../04_polar_vcf
if [ $TARGET = "Rmag" ]; then
    SILWOOD=MAG1,MAG2,MAG3,RMAG402,RMAG403,RMAG404,RMAG407,RMAG707
else
    SILWOOD=AK11,AK15,AK27,RSOC401,RSOC405,RSOC406,RSOC507,RSOC702
fi

# (1)
python map_sites.py $TARGET $MAF $OUT_DIR

# (2)
bcftools view -v snps -m2 -M2 $VCF1 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' - > $OUT_DIR/${TARGET}_snps.tsv
# The expression drops uncalled indel candidates... These sometimes have the same positions as invariant sites.
bcftools query -e 'INFO/INDEL=1' -f '%CHROM\t%POS\t%REF\t%ALT\n' $VCF2 > $OUT_DIR/${COMPARISON}_sites.tsv
python get_anc_states.py \
    $OUT_DIR/${TARGET}_map.tsv \
    $OUT_DIR/${TARGET}_snps.tsv \
    $OUT_DIR/${COMPARISON}_sites.tsv \
    $OUT_DIR/${TARGET}_biSNPs_polarStats.tsv \
    $OUT_DIR/${TARGET}_biSNPs_AA.tsv
bgzip -f $OUT_DIR/${TARGET}_biSNPs_AA.tsv && tabix -f -s1 -b2 -e2 $OUT_DIR/${TARGET}_biSNPs_AA.tsv.gz

# (3)
echo '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral allele">' > $OUT_DIR/aa.hdr
bcftools annotate \
    -a $OUT_DIR/${TARGET}_biSNPs_AA.tsv.gz \
    -c CHROM,POS,INFO/AA \
    -h $OUT_DIR/aa.hdr \
    -Oz -o $OUT_DIR/${TARGET}_polar.vcf.gz \
    $VCF1
bcftools index -f $OUT_DIR/${TARGET}_polar.vcf.gz

# (4) Stats for SNPs polymorphic among the Silwood individuals
bcftools view -s $SILWOOD -v snps -m2 -M2 $VCF1 | bcftools view -i 'AC>0 && AC<AN' | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' - > $OUT_DIR/${TARGET}_snps.silwood.tsv
python get_anc_states.py \
    $OUT_DIR/${TARGET}_map.tsv \
    $OUT_DIR/${TARGET}_snps.silwood.tsv \
    $OUT_DIR/${COMPARISON}_sites.tsv \
    $OUT_DIR/${TARGET}_biSNPs_polarStats.silwood.tsv \
    /dev/null