import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from itertools import combinations
import sys
from sexsigns_utils.calc import get_biallelic_01, to_tskit
from sexsigns_utils.yeast import read_vcf, empty_ts, get_features, intersect_intervals
from sexsigns_utils.sticcs import run_sticcs

import warnings
warnings.filterwarnings("ignore")

##############
### About:
# Run sticcs on phased chromosome alignments.
# Could be modified to alternatively run on unphased data or in the unphased mode (same data; ploidy = 2).
# Run for all clades / chromosomes using run_sticcs_yeast_all.py

##############
### Run as:
# python run_sticcs_yeast.py 1 pair True 12.West_African_cocoa 1
# python run_sticcs_yeast.py 1 pair True '1.Wine_European_(subclade_3)' 4

##############
### Input:
##############
# Command line
ploidy = int(sys.argv[1]) # Controls whether sticcs is run in the phased (ploidy = 1) or unphased (ploidy = 2) mode
strategy = sys.argv[2] # pair; an option to run in 'full' mode (all sample at a time) could be added
sc = sys.argv[3] # True or False; 'second chances' in sticcs
clade = sys.argv[4] # e.g., 12.West_African_cocoa
chrom_id = int(sys.argv[5]) # e.g., 1
##############
# Config
vcf = f"../../0.calls/PH/0.4.total_vcfs/chrom_{chrom_id}/chrom_{chrom_id}.biSNPs.polar.vcf.gz"
callable_df = f"../../config/callable_regions/{clade}.tsv"
callable_df = pd.read_csv(callable_df, sep="\t")
features_df = "../../config/chrom_features.tsv"
features_df = pd.read_csv(features_df, sep="\t")
chrom_len = get_features(features_df, chrom_id)[2]
##############
### Output:
##############
out_dir = Path(f"./ploidy~{ploidy}_strategy~{strategy}_sc~{sc}/chrom_{chrom_id}/{clade}")
out_dir.mkdir(exist_ok=True, parents=True)

##############
### Functions:
def get_callable_regions(callable_df, sample):
    callable_df = callable_df[callable_df["sample"] == sample]
    callable_regions = callable_df[["start", "end"]].to_numpy(dtype=int) 
    return callable_regions

# This can be done more simply
def get_all_regions(callable_regions, chrom_len):
    if len(callable_regions) == 0:
        return [(0, chrom_len, "empty")]
    regions = []
    previous_end = 0
    for start, end in callable_regions:
        if start > previous_end: # Protects from start = 0, but this shouldn't happen
            regions.append((previous_end, start, "empty"))
        regions.append((start, end, "callable"))
        previous_end = end
    if previous_end < chrom_len:
        regions.append((previous_end, chrom_len, "empty"))
    return regions

# Read in data from VCF
genos, pos, samples = read_vcf(vcf)

# Get idx of target samples
targets = callable_df["sample"].unique()
target_idx = np.where(np.isin(samples, targets))[0]

##################
### Infer sticcs ARGs!
# Iterate over pairs of target samples to infer pairwise ARGs
if strategy == "pair":
    for i, j in tqdm(list(combinations(target_idx, r=2))):
    # for i, j in list(itertools.combinations(target_idx, r=2)):
        sample_i = samples[i]
        sample_j = samples[j]

        # Get per-pair callable regions
        callable_regions_i = get_callable_regions(callable_df, chrom_id, sample_i)
        callable_regions_j = get_callable_regions(callable_df, chrom_id, sample_j)
        callable_regions = intersect_intervals(callable_regions_i, callable_regions_j)
        # print(callable_regions)
        # break

        # Get genotypes
        genos_ij = genos[:, [i,j]]
        genos_ij, pos_ij = get_biallelic_01(to_tskit(genos_ij), pos) # Returns tskit-style genos

        ts_lst = []
        regions = get_all_regions(callable_regions, chrom_len)
        for s, e, kind in regions:
            
            if kind == "empty":
                ts = empty_ts(sequence_length=e-s)
                ts_lst.append(ts)
        
            elif kind == "callable":
                mask = (pos_ij >= s) & (pos_ij < e)
                genos_ijk = genos_ij[mask]
                # Normalise positions for the given region
                # (required by sticcs)
                pos_ijk = pos_ij[mask] - s

                ### Validate inputs ###
                if len(genos_ijk) == 0:
                    # print(f"[{clade}; {}~{}; chrom_{chrom_id}; {clade}] Skipping arm {idx+1} in {pair}; no biallelic SNPs!")
                    ts = empty_ts(sequence_length=e-s)
                    ts_lst.append(ts)
                    continue
                # Only regions with non-singleton SNPs
                if np.sum(np.sum(genos_ijk, axis=1) > 1) == 0:
                    # print(f"[chrom_{chrom_id}; {clade}] Skipping arm {idx+1} in {pair}; singletons only!")
                    ts = empty_ts(sequence_length=e-s)
                    ts_lst.append(ts)
                    continue

                # Run sticcs
                ts = run_sticcs(
                    genos_ijk,
                    pos_ijk,
                    e - s,
                    forced_ploidy=ploidy,
                    second_chances=(sc == "True"),
                    silent=True
                )
                ts_lst.append(ts)
        
        # Concatenate into a final ts
        ts = ts_lst[0]
        for ts_k in ts_lst[1:]:
            ts = ts.concatenate(ts_k)
        # print(np.sum([t.span for t in ts.trees() if t.num_roots == 1]) / ts.sequence_length)
        ts.dump(out_dir / f"{sample_i}~{sample_j}.trees")
# elif strategy == "full":