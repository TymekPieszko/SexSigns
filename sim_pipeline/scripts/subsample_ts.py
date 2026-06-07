import tskit, pyslim, msprime, gzip
from pathlib import Path
import numpy as np
from functions_smk import subsample_ts

##########################
### Get paths / params ###
##########################
# Input
ts_in = str(snakemake.input[0])
# Output
ts_out = str(snakemake.output.sub_ts)
vcf_out = str(snakemake.output.vcf)
fasta_dir_out = str(snakemake.output.fasta_dir)
# Other
rep = int(snakemake.wildcards.rep)
mut = float(snakemake.wildcards.mut)
inds_file = str(snakemake.params.inds_file)
L = int(snakemake.params.L)
N_recap = int(snakemake.params.N_recap)
rec_recap = float(snakemake.params.rec_recap)
interval_len = int(snakemake.params.interval_len)
intervals = np.arange(0, L + interval_len, interval_len)
##########################

##########################
### (1) A ts for 2 target inds.
targets = np.loadtxt(inds_file, dtype=int)[rep]
ts = tskit.load(ts_in)
ts = pyslim.recapitate(ts, ancestral_Ne=N_recap, recombination_rate=rec_recap)
ts = msprime.sim_mutations(ts, rate=mut)
ts_2 = subsample_ts(ts, targets)
ts_2.dump(ts_out)

##########################
### (2) A VCO for 2 target inds.
with open(vcf_out, "w") as f:
    ts_2.write_vcf(f, allow_position_zero=True)
f.close()

##########################
### (3) FASTAs for 100 inds (targets included).
Path(fasta_dir_out).mkdir(parents=True, exist_ok=True) # Not sure why Snakemake does not create it!?
rng = np.random.default_rng(rep) # Ensures reproducibility!
inds_100 = np.concatenate([
        targets,
        rng.choice(np.setdiff1d(np.arange(1000), targets), 98, replace=False),
    ]
)
# Annoying but necessary:
# After subsampling, individuals are internally reordered.
# Record the indices of the two target individuals in the FASTA output.
idx_1, idx_2 = np.argsort(np.argsort(inds_100))[0:2]
ts_100 = subsample_ts(ts, inds_100)
for s, e in zip(intervals[:-1], intervals[1:]):
    ts_100_i = ts_100.keep_intervals([[s, e]]).trim()
    with gzip.open(fasta_dir_out + f"/inds-{idx_1}-{idx_2}_span-{s}-{e}.fa.gz", "wt") as f:
        ts_100_i.write_fasta(f)