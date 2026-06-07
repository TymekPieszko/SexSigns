from pathlib import Path
import json
import sys
from sexsigns_utils.calc import calc_category_props
import matplotlib.pyplot as plt

#############
### ABOUT ###
#############
# Plot tree composition under obligate asexuality.

### Run as:
# python plot_tcomp_asex_vPaper.py GC_sim_b~0.1 GC_sim_b~1.0

in_name_1 = sys.argv[1] # Target
in_name_2 = sys.argv[2] # Ground truth to compare against

###################
### Input:
with open(f"../stats/tcomp/{in_name_1}.json", "r") as f:
    data_1 = json.load(f)
with open(f"../stats/tcomp/{in_name_2}.json", "r") as f:
    data_2 = json.load(f)

###################
### Output:
plot_file = Path(f"./tcomp/obligate_asex/{in_name_1}_vs_{in_name_2}.png")
plot_file.parent.mkdir(exist_ok=True, parents=True)

###################
### Get a sorted rec_lst
rec_lst = sorted({float(k.split("_")[1]) for k in data_1.keys()})
print(rec_lst)

###################
### Plot
colors = ["#377eb8", "#4daf4a", "#ff7f00"]
rec_labels = {1e-6:"5.0e-03", 3.162e-07:"1.6e-03", 1e-7:"5.0e-04", 3.162e-08:"1.6e-04", 1e-8:"5.0e-05", 3.162e-09:"1.6e-05", 1e-9:"5.0e-06", 0.0:"0.0"}
fig, axes = plt.subplots(nrows=1, ncols=len(rec_lst), sharey=True, figsize=(40, 8))

for i, rec in enumerate(rec_lst):
    ax = axes[i]
    key = f"0.0_{rec}"

    # Plot target bars
    clarsx_1 = calc_category_props(data_1[key])
    ax.bar(["CL", "AR", "SX"], clarsx_1.values(), width=0.8, color=colors)

    # Plot ground-truth lines
    clarsx_2 = calc_category_props(data_2[key])
    idx = 0
    for tree_cat, prop, color in zip(clarsx_2.keys(), clarsx_2.values(), colors):
        # print(tree_cat, prop, color)
        if tree_cat not in ["AR"]:
            idx += 1
            continue
        ax.hlines(y=prop, xmin=idx-0.4, xmax=idx+0.4, color=color, linestyle=":", linewidth=8)
        idx += 1

    ax.set_ylim(0, 1)
    ax.set_yticks(ticks=[0,0.2,0.4,0.6,0.8,1.0], labels=[0,0.2,0.4,0.6,0.8,1.0], fontsize=32)
    ax.tick_params(axis="x", labelsize=28)
    ax.set_xlabel(rec_labels[rec], fontsize=32, rotation=45)

plt.suptitle(f"{in_name_1} (target) vs {in_name_2} (ground truth)")
plt.tight_layout()
plt.savefig(plot_file)
plt.close()