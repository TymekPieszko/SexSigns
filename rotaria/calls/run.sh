#!/usr/bin/env bash
set -euo pipefail

source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake

# Bioinfo tools (fastqc, fastp, ...) live in the rotspec env
export PATH=/data/biol-bdelloids/scro4331/envs/rotspec/bin:$PATH

snakemake \
    --executor slurm \
    --default-resources slurm_account=biol-bdelloids slurm_partition=short \
    --jobs ${1:-25}
