import tskit, pyslim, msprime, gzip, json
from pathlib import Path
import sys
import numpy as np
from functions_smk import subsample_ts

##########################
### Get paths / params ###
##########################
# Input
ts_in = str(snakemake.input[0])
# Output
ts_out = str(snakemake.output.sub_ts)
# This distinguishes between the 'subsample' and 'subsample_for_abc' instances;
# VCF and FASTA outputs are NOT produced in the latter case (see below). 
with_vcf_fasta = ("vcf" in snakemake.output.keys()) and ("fasta_dir" in snakemake.output.keys())
if with_vcf_fasta:
    vcf_out = str(snakemake.output.vcf)
    fasta_dir_out = str(snakemake.output.fasta_dir)
# Other
rep = int(snakemake.wildcards.rep)
inds_file = str(snakemake.params.inds_file)
if "mut_file" in snakemake.params.keys():
    with open(str(snakemake.params.mut_file)) as f:
        mut = float(json.load(f)[f"{snakemake.wildcards.sex}_{snakemake.wildcards.rec}"])
else:
    mut = float(snakemake.wildcards.mut)
L = int(snakemake.params.L)
N_recap = int(snakemake.params.N_recap)
rec_recap = float(snakemake.params.rec_recap)
interval_len = int(snakemake.params.interval_len)
intervals = np.arange(0, L + interval_len, interval_len)
##########################

##########################
### (1) A ts for the target inds.
targets = np.loadtxt(inds_file, dtype=int, ndmin=2)[rep]
ts = tskit.load(ts_in)
ts = pyslim.recapitate(ts, ancestral_Ne=N_recap, recombination_rate=rec_recap, random_seed=rep + 1)
ts = msprime.sim_mutations(ts, rate=mut, random_seed=rep + 1)
ts_sub = subsample_ts(ts, targets)
ts_sub.dump(ts_out)

### IMPORTANT ###
# If VCF and FASTA dir are NOT in the required outputs,
# exit the script at this point.
if not with_vcf_fasta:
    sys.exit(0)

##########################
### (2) A VCF for 2 target inds.
with open(vcf_out, "w") as f:
    ts_sub.write_vcf(f, allow_position_zero=True)
f.close()

##########################
### (3) FASTAs for 100 inds (targets included).
# Decided to keep it here, but the FASTAs are not used downstream in this version;
# could be used by IQ-TREE, for example.
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
