from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
import tskit, itertools, sys, warnings
from sexsigns_utils.calc import to_tskit, to_allel, get_biallelic_01_allel, get_tree_intervals, calc_delta, count_sgl, calc_kappa, calc_hi, calc_fis
from sexsigns_utils.yeast import read_vcf, get_targets, get_features, get_callable_regions, get_callable_mask, get_callable_length, empty_ts

warnings.simplefilter("ignore")

##############
### Run as:
# python calc_stats.py 12.West_African_cocoa
# python calc_stats.py '1.Wine_European_(subclade_3)'

##############
### Input:
##############
# Command line
clade = sys.argv[1]
clade_dict = "../../../config/clade_dict.json"
# Get target samples
targets = get_targets(clade_dict, clade)
##############
# Config
# VCF defined in the loop
callable_df = f"../../../config/callable_regions/{clade}.tsv"
callable_df = pd.read_csv(callable_df, sep="\t")
features_df = "../../../config/chrom_features.tsv"
features_df = pd.read_csv(features_df, sep="\t")
chrom_ids = list(range(1,17))
chrom_lengths = np.array([get_features(features_df, i)[2] for i in chrom_ids])
# print(chrom_lengths)

###############
### Output:
total_file = Path(f"../stats/{clade}.total.tsv")
chrom_file = Path(f"../stats/{clade}.chrom.tsv")
pair_file = Path(f"../stats/{clade}.pair.tsv")
total_file.parent.mkdir(parents=True, exist_ok=True)
chrom_file.parent.mkdir(parents=True, exist_ok=True)
pair_file.parent.mkdir(parents=True, exist_ok=True)
total_file.write_text("")
chrom_file.write_text("")
pair_file.write_text("")

##############
### Main:

# Initiate summary df
chrom_df = pd.DataFrame(
    index=chrom_ids,
    columns=["delta", "sgl_per_CL", "kappa", "dbl", "trp", "hi", "fis"],
    dtype=float
)

###########################
### Global accumulator matrices
# Accumulate values across chromosomes
# for obtaining global weighted means.
# num = numerator; den = denominator
n = len(targets)
delta_num = np.zeros((n, n), dtype=float)    # sum(delta * span)
delta_den = np.zeros((n, n), dtype=float)    # sum(span)
kappa_num = np.zeros((n, n), dtype=float)    # sum(kappa * chrom_length)
kappa_den = np.zeros((n, n), dtype=float)    # sum(chrom_length)

# For global normalisation of delta (kappa/hi/fis normalised by chrom_length)
CL_spans = []

with open(total_file, "w") as f:
    for chrom_id in tqdm(chrom_ids):
        # Read in data from VCF
        vcf = f"../../../0.calls/PH/0.4.total_vcfs/chrom_{chrom_id}/chrom_{chrom_id}.biSNPs.polar.vcf.gz"

        #######################
        ### Read in data
        genos, pos, samples = read_vcf(vcf)
        chrom_length = chrom_lengths[chrom_id-1]

        #######################
        ### Subset genos to targets
        target_idx = np.where(np.isin(samples, targets))[0] # Get the indices of target samples in this VCF
        target_map = {orig_idx: sub_idx for sub_idx, orig_idx in enumerate(target_idx)}
        ordered_targets = samples[target_idx] # Matches the column order in genos_target
        genos_target = genos[:, target_idx]
        genos_target, pos_target = get_biallelic_01_allel(genos_target, pos)

        #######################
        ### Mask genos
        callable_regions_per_sample = [get_callable_regions(callable_df, chrom_id, target) for target in ordered_targets] # Using samples[target_idx] ensures same order as in genos_target  
        callable_mask_per_sample = get_callable_mask(callable_regions_per_sample, pos_target)
        genos_target.mask = callable_mask_per_sample

        #######################
        ### Calc hi & fis; one value per set of targets
        callable_length = get_callable_length(callable_regions_per_sample)
        hi = calc_hi(genos_target, callable_length)
        fis = calc_fis(genos_target)
        # print(hi, fis)

        ts_dir = f"../../../1.args/sticcs/ploidy~1_strategy~pair_sc~True/chrom_{chrom_id}/{clade}"

        delta_df = pd.DataFrame(index=ordered_targets, columns=ordered_targets, dtype=float)
        sgl_df = pd.DataFrame(index=ordered_targets, columns=ordered_targets, dtype=float)
        kappa_df = pd.DataFrame(index=ordered_targets, columns=ordered_targets, dtype=float)
        dbl_df = pd.DataFrame(index=ordered_targets, columns=ordered_targets, dtype=float)
        trp_df = pd.DataFrame(index=ordered_targets, columns=ordered_targets, dtype=float)

        # For normalising delta on this chromosome
        delta_chrom = 0.0
        CL_chrom = 0.0
        # Iterate pairs of individuals;
        # will be subsetting from genos_target
        for i, j in itertools.combinations(target_idx, 2):
            i = target_map[i]
            j = target_map[j]

            sample_i = ordered_targets[i]
            sample_j = ordered_targets[j]
            pair = f"{sample_i}~{sample_j}"

            genos_ij, pos_ij = get_biallelic_01_allel(genos_target[:, [i, j]], pos_target)

            # Validate input
            if genos_ij.shape[0] == 0:
                continue

            # Calc delta, sgl
            ts_path = Path(ts_dir) / f"{pair}.trees"
            try:
                ts = tskit.load(ts_path)
            except:
                ts = empty_ts(chrom_length)
            interval_dict = get_tree_intervals(ts, include_SX=False) # Get SX intervals
            delta, span = calc_delta(to_tskit(genos_ij), pos_ij, interval_dict)
            sgl = count_sgl(to_tskit(genos_ij), pos_ij, interval_dict) # Mean per interval 
            # For normalising
            if not np.isnan(delta):
                delta_chrom += delta * span
                CL_chrom += span

            # Calc kappa
            genos_ij_tskit = to_tskit(genos_ij)
            genos_ij_tskit = genos_ij_tskit[np.all(genos_ij_tskit != -1, axis=1)] # Exclude sites with missingness
            kappa = calc_kappa(genos_ij_tskit)
            # Below is repeated from calc_kappa():
            der_counts = genos_ij_tskit.sum(axis=1)
            dbl = np.sum(der_counts == 2)
            trp = np.sum(der_counts == 3)

            # Update per-chrom dfs
            delta_df.loc[sample_i, sample_j] = delta
            delta_df.loc[sample_j, sample_i] = delta
            sgl_df.loc[sample_i, sample_j] = sgl
            sgl_df.loc[sample_j, sample_i] = sgl
            kappa_df.loc[sample_i, sample_j] = kappa
            kappa_df.loc[sample_j, sample_i] = kappa
            dbl_df.loc[sample_i, sample_j] = dbl
            dbl_df.loc[sample_j, sample_i] = dbl
            trp_df.loc[sample_i, sample_j] = trp
            trp_df.loc[sample_j, sample_i] = trp

            # Update global accumulators
            if not np.isnan(delta) and span > 0:
                delta_num[i, j] += delta * span
                delta_num[j, i] += delta * span
                delta_den[i, j] += span
                delta_den[j, i] += span
            if not np.isnan(kappa):
                kappa_num[i, j] += kappa * chrom_length
                kappa_num[j, i] += kappa * chrom_length
                kappa_den[i, j] += chrom_length
                kappa_den[j, i] += chrom_length

        # Calculate & record per-chromosome means
        delta_vals = delta_df.values[np.triu_indices_from(delta_df.values, k=1)]
        delta_vals = delta_vals[~np.isnan(delta_vals)]

        sgl_vals = sgl_df.values[np.triu_indices_from(sgl_df.values, k=1)]
        sgl_vals = sgl_vals[~np.isnan(sgl_vals)]

        kappa_vals = kappa_df.values[np.triu_indices_from(kappa_df.values, k=1)]
        kappa_vals = kappa_vals[~np.isnan(kappa_vals)]

        dbl_vals = dbl_df.values[np.triu_indices_from(dbl_df.values, k=1)]
        dbl_vals = dbl_vals[~np.isnan(dbl_vals)]

        trp_vals = trp_df.values[np.triu_indices_from(trp_df.values, k=1)]
        trp_vals = trp_vals[~np.isnan(trp_vals)]

        delta_mean = delta_chrom / CL_chrom if CL_chrom > 0 else np.nan
        sgl_mean = np.nanmean(sgl_vals)
        kappa_mean = np.nanmean(kappa_vals)
        dbl_mean = np.nanmean(dbl_vals)
        trp_mean = np.nanmean(trp_vals)

        chrom_df.loc[chrom_id] = [
            delta_mean,
            sgl_mean,
            kappa_mean,
            dbl_mean,
            trp_mean,
            hi,
            fis,
        ]

        # Collect spans of CL intervals for normalising delta
        CL_spans.append(CL_chrom)

        #####################
        ### Write the per-chromosome pairwise tables into total_file
        col0 = pd.DataFrame({"": [""] * len(delta_df)})
        combined_df = pd.concat(
            [delta_df.reset_index(), col0, sgl_df.reset_index(), col0, kappa_df.reset_index(), col0, dbl_df.reset_index(), col0, trp_df.reset_index()],
            axis=1)
        f.write(f"~ Chromosome {chrom_id} ~\n")
        sep =  "\t" * (delta_df.shape[1] + 2)
        f.write("delta" + sep + "sgl_per_CL" + sep + "kappa" + sep + "dbl" + sep + "trp\n")
        combined_df.to_csv(
            f,
            sep="\t",
            float_format="%.6g",
            index=False)
        f.write("\n")
        #####################


delta_means = chrom_df["delta"].values.astype(float)
kappa_vals = chrom_df["kappa"].values.astype(float)
hi_vals = chrom_df["hi"].values.astype(float)
fis_vals = chrom_df["fis"].values.astype(float)

#######################
### Delta mean across chromosomes (weight by CL span)
CL_spans = np.array(CL_spans, dtype=float)
delta_mask = ~np.isnan(delta_means)
delta_mean = (
    np.sum(delta_means[delta_mask] * CL_spans[delta_mask]) / np.sum(CL_spans[delta_mask])
    if np.any(delta_mask) else np.nan
)
#######################

#######################
### Kappa/hi/fis means (weight by chromosome length)
chrom_lengths = np.array(chrom_lengths, dtype=float)

kappa_mask = ~np.isnan(kappa_vals)
kappa_mean = (
    np.sum(kappa_vals[kappa_mask] * chrom_lengths[kappa_mask]) / np.sum(chrom_lengths[kappa_mask])
    if np.any(kappa_mask) else np.nan
)

hi_mask = ~np.isnan(hi_vals)
hi_mean = (
    np.sum(hi_vals[hi_mask] * chrom_lengths[hi_mask]) / np.sum(chrom_lengths[hi_mask])
    if np.any(hi_mask) else np.nan
)

fis_mask = ~np.isnan(fis_vals)
fis_mean = (
    np.sum(fis_vals[fis_mask] * chrom_lengths[fis_mask]) / np.sum(chrom_lengths[fis_mask])
    if np.any(fis_mask) else np.nan
)
#######################

chrom_df.loc["mean"] = [
    delta_mean,
    np.nan,
    kappa_mean,
    np.nan,
    np.nan,
    hi_mean,
    fis_mean,
]

# Compute per-pair means across chromosomes
delta_pair_means = np.where(delta_den > 0, delta_num / delta_den, np.nan)
kappa_pair_means = np.where(kappa_den > 0, kappa_num / kappa_den, np.nan)

# Long table
a_idx, b_idx = np.triu_indices(n, k=1)
targets_arr = np.array(ordered_targets)
pair_df = pd.DataFrame({
    "ind_A": targets_arr[a_idx],
    "ind_B": targets_arr[b_idx],
    "delta": delta_pair_means[a_idx, b_idx],
    "kappa": kappa_pair_means[a_idx, b_idx],
}).set_index(["ind_A", "ind_B"])

############################
### Write summary files!
with open(chrom_file, "a") as f:
    f.write("# Note:\n")
    f.write("# delta is normalised by total CL span per pair, then by total CL span per chromosome\n")
    f.write("# kappa / hi / fis are normalised by chromosome lengths\n")
    f.write("# sgl_per_CL, dbl, trp are NOT normalised (mean across pairs)\n\n")
    f.write("# Per-chromosome summary (rows = chromosomes; mean = global normalised mean)\n")
    chrom_df.to_csv(f, sep="\t", float_format="%.6g")
with open(pair_file, "a") as f:
    f.write("# Per-pair summary (across chromosomes)\n")
    pair_df.to_csv(f, sep="\t", float_format="%.6g")