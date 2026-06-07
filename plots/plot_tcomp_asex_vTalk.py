from pathlib import Path
import json
import sys
from sexsigns_utils.calc import calc_category_props
import matplotlib.pyplot as plt

#############
### ABOUT ###
#############
# Plot tree composition under obligate asexuality.

# python plot_tcomp_asex_vTalk.py GC_sim_b~1.0

in_name = sys.argv[1]

###################
### Input:
with open(f"../stats/tcomp/{in_name}.json", "r") as f:
    data = json.load(f)

###################
### Output:
plot_file = Path(f"./tcomp/obligate_asex/{in_name}.png")
plot_file.parent.mkdir(exist_ok=True, parents=True)

###################
### Get a sorted rec_lst
rec_lst = sorted({float(k.split("_")[1]) for k in data.keys()})
print(rec_lst)

###################
### Plot
colors = ["#377eb8", "#4daf4a", "#ff7f00"]
rec_labels = {1e-6:"5.0e-03", 3.162e-07:"1.6e-03", 1e-7:"5.0e-04", 3.162e-08:"1.6e-04", 1e-8:"5.0e-05", 3.162e-09:"1.6e-05", 1e-9:"5.0e-06", 0.0:"0.0"}
fig, axes = plt.subplots(nrows=1, ncols=len(rec_lst), sharey=True, figsize=(40, 8))

for i, rec in enumerate(rec_lst):
    ax = axes[i]
    key = f"0.0_{rec}"

    clarsx_1 = calc_category_props(data[key])
    bars = ax.bar(["CL", "AR", "SX"], clarsx_1.values(), width=0.8, color=colors)
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,
            height,
            f"{height:.3f}",
            ha="center",
            va="bottom",
            fontsize=24
        )
    
    # Y axis
    ax.set_ylim(0, 1)
    ax.set_yticks(ticks=[0,0.2,0.4,0.6,0.8,1.0], labels=[0,0.2,0.4,0.6,0.8,1.0], fontsize=32)

    # X axis
    ax.tick_params(axis="x", labelsize=28)
    ax.set_xlabel(rec_labels[rec], fontsize=32, rotation=45)

plt.suptitle(in_name)
plt.tight_layout()
plt.savefig(plot_file)
plt.close()