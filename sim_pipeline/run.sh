# Reset trick:
# rm -R .snakemake/ && snakemake --touch

# Create envs:
# module load Mamba/4.14.0-0
# snakemake --use-conda --conda-frontend mamba --conda-create-envs-only
# Might be required: conda config --set channel_priority flexible

# Dry run:
# snakemake -n --cores 1 --rerun-triggers mtime
# snakemake -n --cores 1 --rerun-triggers mtime --forcerun subsample

### Real run
snakemake \
--executor slurm \
-j 75 \
--resources mem_mb=24000 \
--default-resources slurm_account=biol-bdelloids slurm_partition=devel \
--group-components subsample_for_abc=40 \
--use-conda --conda-frontend mamba \
--rerun-triggers mtime \
--rerun-incomplete \
--latency-wait 10000 \
--keep-going
# --forcerun subsample

# --set-resources IQ_TREE:slurm_partition=short 