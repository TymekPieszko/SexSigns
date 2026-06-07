from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
import tskit, tqdm, sys

#############
### ABOUT ###
#############
# Calculate stats for tree sequences subsampled in the simulation pipeline.

### Run as:
# python calc_sub_ts.py GC

model = sys.argv[1]

out_dir = Path("./pipeline_stats/")
out_dir.mkdir(exist_ok=True)

total_blocks = defaultdict(dict)
total_SNVs = defaultdict(dict)

for dir in Path(f"../sim_pipeline_2/sim_output/{model}/1.0.sub_ts").glob(
    "SEX~*/REC~*/MUT~*/"
):

    SEX = float(dir.parts[-3].split("~")[1])
    REC = float(dir.parts[-2].split("~")[1])

    print(SEX, REC)

    blocks_ts = []
    SNVs_ts = []

    for ts in tqdm.tqdm(dir.glob("*trees")):

        ts = tskit.load(ts)
        # print(ts.num_trees, ts.num_sites)

        blocks_ts.append(ts.num_trees)
        SNVs_ts.append(ts.num_sites)

    print(len(blocks_ts), len(SNVs_ts))

    total_blocks[SEX][REC] = np.nanmean(blocks_ts)
    total_SNVs[SEX][REC] = np.nanmean(SNVs_ts)

###################
# Convert to dataframes

total_blocks = (
    pd.DataFrame(total_blocks)
    .sort_index()
    .sort_index(axis=1)
    .iloc[:, :-1][::-1]
)

total_SNVs = (
    pd.DataFrame(total_SNVs)
    .sort_index()
    .sort_index(axis=1)
    .iloc[:, :-1][::-1]
)

sex_labels_in_N_tsv = [
    "0",
    "0.01 N^-1",
    "0.03162 N^-1",
    "0.1 N^-1",
    "0.3162 N^-1",
    "1 N^-1",
    "3.162 N^-1",
    "10 N^-1",
    "31.62 N^-1",
    "100 N^-1",
    "316.2 N^-1",
]

loh_labels_in_N_tsv = [
    "0",
    "0.005 N^-1",
    "0.0158 N^-1",
    "0.05 N^-1",
    "0.158 N^-1",
    "0.5 N^-1",
    "1.58 N^-1",
    "5 N^-1",
]

# Replace numeric labels with text labels

sex_to_label = sorted(total_blocks.columns)
rec_to_label = sorted(total_blocks.index)[::-1]

sex_label_dict = {
    sex: label
    for sex, label in zip(sex_to_label, sex_labels_in_N_tsv)
}

rec_label_dict = {
    rec: label
    for rec, label in zip(rec_to_label, loh_labels_in_N_tsv[::-1])
}

total_blocks.columns = [
    sex_label_dict[v]
    for v in total_blocks.columns
]

total_blocks.index = [
    rec_label_dict[v]
    for v in total_blocks.index
]

total_SNVs.columns = [
    sex_label_dict[v]
    for v in total_SNVs.columns
]

total_SNVs.index = [
    rec_label_dict[v]
    for v in total_SNVs.index
]

###################

total_blocks.to_csv(
    out_dir + f"{model}_blocks.txt",
    sep="\t",
    header=True,
    index=True
)

total_SNVs.to_csv(
    out_dir + f"{model}_SNVs.txt",
    sep="\t",
    header=True,
    index=True
)