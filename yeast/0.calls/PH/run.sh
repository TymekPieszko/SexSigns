# module load Mamba/4.14.0-0

### If snakemake wants to regenerate outputs, do:
# rm -R .snakemake/ && snakemake --touch

### Create envs
# snakemake --use-conda --conda-frontend mamba --conda-create-envs-only


### Dry run
# snakemake -n --cores 1 --rerun-triggers mtime

### Real run
snakemake \
--executor slurm \
-j 80 \
--default-resources slurm_account=biol-bdelloids slurm_partition=devel \
--use-conda --conda-frontend mamba \
--rerun-triggers mtime \
--rerun-incomplete \
--latency-wait 60 \
# --slurm-efficiency-report \
# --forceall

# --resources threads=16 mem_mb=$((5000 * 2)) runtime=10 \
# --set-resources IQ_TREE:slurm_partition=short SINGER:slurm_partition=short \
# --group-components subsample=40 iqtree=120 rootprune=80 singer=80 \
# --forcerun root_prune