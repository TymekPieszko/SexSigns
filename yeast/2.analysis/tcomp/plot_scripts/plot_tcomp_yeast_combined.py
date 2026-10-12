from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json, sys, tqdm
from collections import Counter
from sexsigns_utils.calc import codes, calc_category_props
from sexsigns_utils.yeast import plot_grid

##############
### About:
# Plot tree composition
# Tymek Pieszko; 2026-04-14

##############
### Run as:
# python plot_tcomp_yeast_combined.py

###############
### Input:
clade_dict = f"../../../config/clade_dict.json"
life_cycle_dict = "../../../config/life_cycle_dict.json"
data_dir = Path(f"../data/ploidy~1_strategy~pair_sc~True/")
subsample_brazilian_bioethanol = False

###############
### Output:
plot_dir = Path(f"../plots/ploidy~1_strategy~pair_sc~True/combined")
plot_dir.mkdir(exist_ok=True, parents=True)

###############
### Load data:
with open(life_cycle_dict, 'r') as f:
    life_cycle_dict = json.load(f)
###############
with open(clade_dict, 'r') as f:
    clade_dict = json.load(f)
clades = list(clade_dict.keys())
###############

###############
### Plot 1
cats = ["CL", "AR", "SX"]
colors = ["#377eb8", "#4daf4a", "#ff7f00"]
for clade in tqdm.tqdm(clades):
    with open(data_dir / f"chrom_1/{clade}.json") as f:
        pairs = list(json.load(f).keys())
    total = {pair: {code: 0.0 for code in codes} for pair in pairs}
    for tcomp_per_pair in data_dir.rglob(f"*/{clade}.json"):
        with open(tcomp_per_pair, 'r') as f:
            tcomp_per_pair = json.load(f)
        for pair, tcomp in tcomp_per_pair.items():
            for code, span in tcomp.items():
                total[pair][code] += span
    props_arr = np.array([[calc_category_props(t)[c] for c in cats] for t in total.values()], dtype=float)
    mins, maxs = props_arr.min(axis=0), props_arr.max(axis=0)
    pooled = {code: sum(t[code] for t in total.values()) for code in codes}
    pooled_props = np.array([calc_category_props(pooled)[c] for c in cats])
    fig = plt.figure(figsize=(6.4, 6.9))
    ax = fig.add_axes([1.5 / 6.4, 0.4 / 6.9, 4.5 / 6.4, 6 / 6.9])
    for i in range(len(cats)):
        x = [i, i + 1]
        ax.fill_between(x=x, y1=mins[i], y2=maxs[i], color=colors[i], alpha=0.2, linewidth=0)
        ax.plot(x, [mins[i], mins[i]], color=colors[i], linewidth=6)
        ax.plot(x, [maxs[i], maxs[i]], color=colors[i], linewidth=6)
        ax.plot(x, [pooled_props[i], pooled_props[i]], color=colors[i], linestyle=(0, (1, 0.5)), linewidth=14)
    ax.set_xticks(np.arange(len(cats)) + 0.5)
    ax.tick_params(axis="x", labelbottom=False, length=18, width=3)
    ax.tick_params(axis="y", labelsize=38)
    ax.set_xlim(0, len(cats))
    ax.set_ylim(0, 1)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    plot_file = plot_dir / f"{clade}.total.png"
    fig.savefig(plot_file, dpi=150)
    plt.close()
###############
### Plot 2
# Get pairs
for clade in tqdm.tqdm(clades):
    
    # Get pairs of individuals
    with open(data_dir / f"chrom_1/{clade}.json") as f:
        tcomp_per_pair = json.load(f)
    pairs = list(tcomp_per_pair.keys())
    # This is specific to Brazilian Bioethanol
    if subsample_brazilian_bioethanol:
        to_remove = ["CMT", "CNG", "CNK", "CNM", "CNP", "CNR", "BVB", "CNI"]
        pairs = [pair for pair in pairs if not (pair.split("~")[0] in to_remove or pair.split("~")[1] in to_remove)]

    # Get samples
    samples = clade_dict[clade]
    if subsample_brazilian_bioethanol:
        samples = [s for s in samples if s not in to_remove]

    # Combine 
    total = {pair: {code: 0.0 for code in codes} for pair in pairs}
    for tcomp_per_pair in data_dir.rglob(f"*/{clade}.json"):
        # print(tcomp_per_pair)
        with open(tcomp_per_pair, 'r') as f:
            tcomp_per_pair = json.load(f)
        for pair, tcomp in tcomp_per_pair.items():
            if pair not in total:
                continue
            for code, span in tcomp.items():
                total[pair][code] += span
    for pair in total:
        s = sum(total[pair].values())
        if s > 0:
            total[pair] = {k: v / s for k, v in total[pair].items()}
    
    # Plot
    max_cats = []
    means = {code:[] for code in ["CL", "AR", "SX"]}
    for tcomp in total.values():
        clarsx = calc_category_props(tcomp)
        max_code = max(clarsx, key=clarsx.get)
        max_cats.append(max_code)
        means[max_code].append(clarsx[max_code])
    print(clade, Counter(max_cats))
    means = {k: np.mean(v) for k,v in means.items()}
    print(means)
    fig, ax = plot_grid(total, life_cycle_dict, samples, figsize_scale=1.2, font_size=0)
    plot_file = plot_dir / f"{clade}.pairwise.png" 
    # plt.suptitle(clade, fontsize=50, x=0.50, y=0.98, ha="left")
    plt.savefig(plot_file)
    plt.close()