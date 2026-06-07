from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import sys
from sexsigns_utils.yeast import intersect_intervals, get_features

# python calc_callable_props.py 12.West_African_cocoa
# python calc_callable_props.py '1.Wine_European_(subclade_3)'
# python calc_callable_props.py 3.Brazilian_bioethanol
# python calc_callable_props.py 25.Sake

clade = sys.argv[1]

callable_file = Path(f"./callable_regions/{clade}.tsv")
features = pd.read_csv("../config/chrom_features.tsv", sep="\t")
df = pd.read_csv(callable_file, sep="\t")

def interval_length(intervals):
    if len(intervals) == 0:
        return 0
    return int(np.sum(intervals[:, 1] - intervals[:, 0]))

# Total genome length
genome_len = 0
for chrom_id in range(1, 17):
    genome_len += get_features(features, chrom_id)[2]

# Organize intervals by sample and chromosome
samples = sorted(df["sample"].unique())

by_sample = {
    sample: {
        chrom_id: subdf[["start", "end"]].to_numpy(dtype=int)
        for chrom_id, subdf in sdf.groupby("chrom_id")
    }
    for sample, sdf in df.groupby("sample")
}

# -------------------
# Per-pair output
# -------------------
pair_rows = []

for s1, s2 in combinations(samples, 2):
    shared_bases = 0

    for chrom_id in range(1, 17):
        intervals_1 = by_sample.get(s1, {}).get(chrom_id, np.empty((0, 2), dtype=int))
        intervals_2 = by_sample.get(s2, {}).get(chrom_id, np.empty((0, 2), dtype=int))

        shared = intersect_intervals(intervals_1, intervals_2)
        shared_bases += interval_length(shared)

    pair_rows.append({
        "sample_1": s1,
        "sample_2": s2,
        "shared_callable_bases": shared_bases,
        "shared_callable_prop": shared_bases / genome_len,
    })

pair_out = pd.DataFrame(pair_rows)

pair_mean_row = pd.DataFrame([{
    "sample_1": "mean",
    "sample_2": "",
    "shared_callable_bases": pair_out["shared_callable_bases"].mean() if len(pair_out) else np.nan,
    "shared_callable_prop": pair_out["shared_callable_prop"].mean() if len(pair_out) else np.nan,
}])

pair_out = pd.concat([pair_out, pair_mean_row], ignore_index=True)

pair_out_path = Path(f"./callable_props/{clade}.callable_props.per_pair.tsv")
pair_out.to_csv(pair_out_path, sep="\t", index=False)

# -------------------
# Per-sample output
# -------------------
sample_rows = []

for sample in samples:
    callable_bases = 0

    for chrom_id in range(1, 17):
        intervals = by_sample.get(sample, {}).get(
            chrom_id, np.empty((0, 2), dtype=int)
        )
        callable_bases += interval_length(intervals)

    sample_rows.append({
        "sample": sample,
        "callable_bases": callable_bases,
        "callable_prop": callable_bases / genome_len,
    })

sample_out = pd.DataFrame(sample_rows)

sample_mean_row = pd.DataFrame([{
    "sample": "mean",
    "callable_bases": sample_out["callable_bases"].mean() if len(sample_out) else np.nan,
    "callable_prop": sample_out["callable_prop"].mean() if len(sample_out) else np.nan,
}])

sample_out = pd.concat([sample_out, sample_mean_row], ignore_index=True)

sample_out_path = Path(f"./callable_props/{clade}.callable_props.per_sample.tsv")
sample_out.to_csv(sample_out_path, sep="\t", index=False)

print(f"Wrote {pair_out_path}")
print(f"Wrote {sample_out_path}")