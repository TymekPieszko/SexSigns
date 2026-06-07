from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json, sys
from sexsigns_utils.calc import codes
from sexsigns_utils.yeast import plot_grid, plot_category_props, calc_spatial_weights, plot_spatial_weights, color_dict

##############
### About:
# Plot tree composition
# Run for all clades / chromosomes using plot_tcomp_yeast_all.py
# Tymek Pieszko; 2026-04-14

##############
### Run as:
# python plot_tcomp_yeast.py 1 pair True 12.West_African_cocoa 1

###############
### Input:
ploidy = int(sys.argv[1])
strategy = sys.argv[2]
sc = sys.argv[3]
clade = sys.argv[4]
chrom_id = int(sys.argv[5])
# Config
clade_dict = f"../../../config/clade_dict.json"
life_cycle_dict = "../../../config/life_cycle_dict.json"
# Data to plot
data_dir = Path(f"../data/ploidy~{ploidy}_strategy~{strategy}_sc~{sc}/chrom_{chrom_id}/")
tcomp_per_pair = data_dir / f"{clade}.json" # A dictionary of dictionaries; tcomp for each pair of isolates
stcomp_per_pair = data_dir / f"{clade}.spatial.json" # A dictionary of dictionaries; stcomp for each pair of isolates
# Get centromere position position for plotting
features = "../../../config/chrom_features.tsv"
features = pd.read_csv(features, sep="\t")
features = features[features["chrom"] == f"chrom_{chrom_id}"]
centr_start, centr_end = features.iloc[0,:][1:3]
centr_mid = (centr_start + centr_end) / 2

###############
### Output:
plot_dir = Path(f"../plots/ploidy~{ploidy}_strategy~{strategy}_sc~{sc}/chrom_{chrom_id}/")
plot_dir.mkdir(exist_ok=True, parents=True)
plot_1 = plot_dir /f"{clade}.total.png"
plot_2 = plot_dir /f"{clade}.pairwise.png"
plot_3 = plot_dir /f"{clade}.spatial.png"

###############
### Load data:
with open(clade_dict, 'r') as f:
    clade_dict = json.load(f)
samples = clade_dict[clade]
# print(samples)
###############
with open(tcomp_per_pair, 'r') as f:
    tcomp_per_pair = json.load(f)
###############
with open(stcomp_per_pair, 'r') as f:
    stcomp_per_pair = json.load(f)
###############
with open(life_cycle_dict, 'r') as f:
    life_cycle_dict = json.load(f)
###############

###############
### Plot 1
total = {code: 0.0 for code in codes}
for tcomp in tcomp_per_pair.values():
    for pair in tcomp.keys():
        total[pair] += tcomp[pair]
total_sum = sum(total.values())
total = {k : v/total_sum for k, v in total.items()}
fig, ax = plt.subplots(figsize=(9,10))
plot_category_props(ax, total, "sample1", "sample2", font_size=24)
plt.xticks(range(3), ["CL", "AR", "SX"], fontsize=50)
ax.tick_params(axis='x', pad=15)
plt.yticks(fontsize=42)
plt.savefig(plot_1)
plt.close()

###############
### Plot 2
fig, ax = plot_grid(total, life_cycle_dict, samples, figsize_scale=1.2, font_size=6)
plt.suptitle(f"ploidy~{ploidy}_strategy~{strategy}_sc~{sc}", fontsize=42, x=0.50, y=0.98, ha="left")
plt.savefig(plot_2)
plt.close()

###############
### Plot 3
stcomp_per_pair = list(stcomp_per_pair.values())
blocks, weights = calc_spatial_weights(stcomp_per_pair, codes + ["poly", "multi_root"])
# print(next((i for i,(s,e) in enumerate(blocks) if s <= 449711 < e), None))
# print([d["multi_root"] for d in weights][425])
plot_spatial_weights(blocks, weights, color_dict, figsize=(12,4))
# plt.xticks(np.arange(0,170000,20000), labels=np.arange(0,17,2), fontsize=17)
# plt.yticks([0.0,0.5,1.0],[0.0,0.5,1.0],fontsize=17)
# plt.xlim(440000, 460000)
plt.axvline(centr_mid, color="black", linewidth=6)
plt.tight_layout()
plt.savefig(plot_3)
plt.close()