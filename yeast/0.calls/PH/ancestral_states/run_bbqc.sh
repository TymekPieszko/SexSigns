#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=100g
#SBATCH --ntasks-per-node=16
#SBATCH --time=00:10:00
#SBATCH --partition=devel
#SBATCH --job-name=bbqc
#SBATCH -o /data/zool-barralab/scro4331/job_stdout/job_%x_%j.out
set -euo pipefail

module load Anaconda3/2022.05
source activate bbtools

SAMPLE=$1
THREADS=12

IN_READS1=./reads/${SAMPLE}_1.fastq.gz
IN_READS2=./reads/${SAMPLE}_2.fastq.gz
TMP_READS1=./reads/${SAMPLE}_1.tmp.fastq.gz
TMP_READS2=./reads/${SAMPLE}_2.tmp.fastq.gz
OUT_READS1=./reads/${SAMPLE}_1.bbduk.ecc.fastq.gz
OUT_READS2=./reads/${SAMPLE}_2.bbduk.ecc.fastq.gz

# Get kmer dist of input
# kmercountexact.sh -Xmx60g t=$THREADS \
#   in1=$READS1 in2=$READS2 \
#   khist=$OUT_PATH/${PREFIX}.khist

# Run bbduk.sh
bbduk.sh -Xmx60g t=$THREADS \
  in1=$IN_READS1 in2=$IN_READS2 \
  out1=$TMP_READS1 out2=$TMP_READS2 \
  ref=./reads/adapters.fasta \
  ktrim=r k=23 mink=11 hdist=2 maq=10 minlen=100 tpe tbo \
  stats=./reads/${SAMPLE}.bbduk.contaminants

# Run tadpole.sh
tadpole.sh -Xmx60g t=$THREADS \
  mode=correct \
  in1=$TMP_READS1 in2=$TMP_READS2 \
  out1=$OUT_READS1 out2=$OUT_READS2 \
  tossjunk

# Get kmer dist of output
# kmercountexact.sh -Xmx60g t=$THREADS \
#   in1=$OUT_PATH/$ECC1 in2=$OUT_PATH/$ECC2 \
#   khist=$OUT_PATH/${PREFIX}.bbduk.ecc.khist

rm ./reads/*tmp*
