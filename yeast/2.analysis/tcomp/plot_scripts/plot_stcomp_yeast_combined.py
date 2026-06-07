from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json, sys, tqdm
from sexsigns_utils.calc import codes
from sexsigns_utils.yeast import calc_spatial_weights, color_dict

##############
### About:
# -------------------
# Tymek Pieszko; 2026-04-14

##############
### Run as:
# python plot_stcomp_yeast_combined.py 1 pair True 0.0 12.West_African_cocoa
# python plot_stcomp_yeast_combined.py 1 pair True 0.0 3.Brazilian_bioethanol
# python plot_stcomp_yeast_combined.py 1 pair True 0.0 '1.Wine_European_(subclade_3)'
# python plot_stcomp_yeast_combined.py 1 pair True 0.0 25.Sake

###############
### Input:
ploidy = int(sys.argv[1])
strategy = sys.argv[2]
sc = sys.argv[3]
min_pop_prop = float(sys.argv[4]) # At a given chromosomal position; threshold prop of sample pairs for plotting tcomp; otheriwse undefined  
clade = sys.argv[5]
features = "../../../config/chrom_features.tsv"
features = pd.read_csv(features, sep="\t")

###############
### Output:
out_dir = Path(f"../plots/ploidy~{ploidy}_strategy~{strategy}_sc~{sc}/combined/")
out_dir.mkdir(parents=True, exist_ok=True)
plot = out_dir / f"{clade}.spatial.undefined_below_{min_pop_prop}.png"

###############
### Functions:
def plot_stcomp_chrom(offset, blocks, weights, color_dict, min_pop_prop):
    for (s, e), w in zip(blocks, weights):
        s += offset
        e += offset 
        bottom = 0
        if min_pop_prop == 0:
            # Do nothing; plot tcomp as is
            pass
        elif w["poly"] + w["multi_root"] > (1-min_pop_prop):
            # Prop of poly + multi_root exceeds threshold; 
            # tcomp undefined.
            w = {"undefined": 1.0}
        else:
            # Prop of poly + multi_root below threshold; 
            # plot tcomp but remove poly + multi_root
            del w["poly"]
            del w["multi_root"]
            total = sum(w.values()) # Will never be 0
            w = {k : v / total for k,v in w.items()}
        for code, prop in w.items():
            if prop == 0:
                continue
            color = color_dict[code]
            hatch = "xxx" if color == "white" else ""
            plt.bar(
                x=s,
                height=prop,
                width=e - s,
                bottom=bottom,
                align='edge',
                color=color,
                hatch=hatch,
                linewidth=0
            )
            bottom += prop
    return blocks[-1][1]

def read_stcomp_per_pair(stcomp_per_pair):
    with open(stcomp_per_pair, 'r') as f:
        stcomp_per_pair = json.load(f)
    return stcomp_per_pair

###############
### Plot:

plt.figure(figsize=(22, 4))
offset = 0
boundries = []
chrom_lengths = []
for chrom_id in tqdm.tqdm(range(1,17)):
    stcomp_per_pair = f"../data/ploidy~{ploidy}_strategy~{strategy}_sc~{sc}/chrom_{chrom_id}/{clade}.spatial.json" # A dictionary of dictionaries; stcomp for each pair of isolates
    stcomp_per_pair = read_stcomp_per_pair(stcomp_per_pair)
    stcomp_per_pair = list(stcomp_per_pair.values()) # Ignore the pairs themselves; not relevant for plotting.
    blocks, weights = calc_spatial_weights(stcomp_per_pair, codes + ["poly", "multi_root"])
    chrom_length = plot_stcomp_chrom(offset, blocks, weights, color_dict, min_pop_prop)

    # Draw centromere
    row = features[features["chrom"] == f"chrom_{chrom_id}"]
    centr_start, centr_end = row.iloc[0,:][1:3].astype(int)
    centr_mid = (centr_start + centr_end) / 2
    plt.axvline(centr_mid+offset, color="darkgray", linewidth=4)

    # Draw chromosome boundry
    plt.axvline(offset, color="black", linewidth=5, linestyle="--")

    # Update position
    boundries.append(offset)
    chrom_lengths.append(chrom_length)
    offset += chrom_length
boundries.append(offset)

# Label chromosomes
boundries = np.array(boundries)
chrom_mids = (boundries[:-1] + boundries[1:]) / 2
chrom_labels = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
for i, pos in enumerate(chrom_mids):
    plt.text(
        pos,
        1.01,
        chrom_labels[i],
        ha="center",
        va="bottom",
        fontsize=24,
        clip_on=False
    )
    plt.text(
        pos,
        -0.03,                     # relative to axis height
        f"{chrom_lengths[i] / 1e6:.2f}",
        ha="center",
        va="top",
        fontsize=20,
        transform=plt.gca().get_xaxis_transform(),
        clip_on=False
    )

# Finish plot
plt.axvline(offset, color="black", linewidth=5, linestyle="--")
plt.xlim(0, offset)
plt.xticks([])
plt.ylim(0, 1)
plt.yticks(fontsize=20)

plt.tight_layout()
plt.subplots_adjust(top=0.90)

plt.savefig(plot, dpi=300)
plt.close()