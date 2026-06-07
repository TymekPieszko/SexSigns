# Read quality control
sbatch run_bbqc.sh AMH
sbatch run_bbqc.sh CEG
sbatch run_bbqc.sh CEI

# Align to reference;
# run this interactively.
aconda minimap2
REF=/data/biol-bdelloids/scro4331/SexSigns_2025/yeast/0.calls/PH/0.0.data/GCA_000146045.2_R64_genomic.2026-01-30.simpleIDs.fasta
SAMPLE=CEI # AMH / CEG / CEI
READS1=./reads/${SAMPLE}_1.bbduk.ecc.fastq.gz
READS2=./reads/${SAMPLE}_2.bbduk.ecc.fastq.gz
BAM=./aligned/${SAMPLE}.sorted.bam
minimap2 -ax sr -t 10 $REF $READS1 $READS2 | \
samtools sort -@ 6 -o $BAM && \
samtools index $BAM

# Call variants;
# bcftools reheader was only necessary because you did not add read groups during mapping!
# You also needed to rename contigs using bcftools annotate --rename-chrs chrom_names.txt outgroups.vcf.gz -Oz -o outgroups_2.vcf.gz;
# this is because you actually DID NOT use the 'simpleIDs' version of the fasta.
# The 'contigs' are now chromosome1, chromosome2, etc.
module load BCFtools/1.14-GCC-11.2.0
module load SAMtools/1.14-GCC-11.2.0
REF=/data/biol-bdelloids/scro4331/SexSigns_2025/yeast/0.calls/PH/0.0.data/GCA_000146045.2_R64_genomic.2026-01-30.simpleIDs.fasta
BAM1=./aligned/AMH.sorted.bam
BAM2=./aligned/CEG.sorted.bam
BAM3=./aligned/CEI.sorted.bam
VCF=./outgroups.vcf.gz
bcftools mpileup -f $REF --threads 12 -Ou $BAM1 $BAM2 $BAM3 | \
bcftools call -mv --threads 12 -Ou | \
bcftools reheader -s samples.txt | \
bcftools view - -i 'QUAL>=30 && TYPE="snp"' -Oz -o $VCF

# Extract chromosomes
for i in {1..16}; do
    bcftools view -r chromosome${i} ./outgroups.vcf.gz -Oz -o ./outgroup_vcfs/chrom_${i}.vcf.gz
    bcftools index ./outgroup_vcfs/chrom_${i}.vcf.gz
done

# Get consensus sequence
# REF_CHROM1=/data/biol-bdelloids/scro4331/SexSigns_2025/yeast/0.calls/PH/0.0.data/GCA_000146045.2_R64_genomic.2026-01-30.chrom1.fasta
# VCF_CHROM1=./outgroups.chrom1.vcf.gz
# bcftools consensus -f $REF_CHROM1 $VCF_CHROM1 -o chrom1.consensus.fasta
