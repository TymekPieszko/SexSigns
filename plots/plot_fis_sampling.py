from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import sys, json
from sexsigns_utils.calc import calc_akaike_weights

#############
### ABOUT ###
#############
# Plot the distribution of F_IS values under different sampling regimes.

### Run as:
# python plot_fis_sampling.py GC
# python plot_fis_sampling.py CO

###################
### Input:
model = sys.argv[1]
in_data = f"../stats/fis_sampling/{model}.json"
n_lst = [5, 10, 15, 50, 100]
if model == "GC":
    rec_min = 3.162e-08
    rec_max = 3.162e-07
if model == "CO":
    rec_min = 6.326010734420079e-10
    rec_max = 6.344556165159151e-09
colors = ["blue", "black", "red"]
gammas = [r"$0.16N^{-1}$", r"$0.5N^{-1}$", r"$1.6N^{-1}$"]

###################
### Output:
plot_file = f"./fis_sampling/{model}.png"
Path(plot_file).parent.mkdir(exist_ok=True, parents=True)

# Read in data & extract info
in_data = json.load(open(in_data))
keys = [k for k in in_data if k.startswith(f"5_")]
recs = np.array([float(k.split("_")[1]) for k in keys])
recs = np.sort(recs)
# Get rec mask
rec_mask = (recs >= rec_min) & (recs <= rec_max)

fig, axes = plt.subplots(len(n_lst), 1, figsize=(6, 10), sharex=True, sharey=True)

for i, n in enumerate(n_lst):
    ax = axes[i]
    dists = [in_data[f"{n}_{rec}"] for rec in recs]
    weights = calc_akaike_weights(obs=0, dists=dists, n_params=1)
    # Dict to use for plotting
    weights = {rec:weight for rec,weight in zip(recs, weights)} 

    for j, rec in enumerate(recs[rec_mask]):
        vals = in_data[f"{n}_{rec}"]
        ax.hist(vals, bins=20, density=True, histtype="step", color=colors[j])
    
    # Mark fis = 0
    ax.axvline(0, linestyle="--", color="black")

    ax.text(0.02, 0.95, rf"$n = {n}$", transform=ax.transAxes, ha="left", va="top", fontsize=16)
    ax.set_yticks([])


    # Add Akaike weights
    for j, (gamma, rec) in enumerate(zip(gammas, recs[rec_mask])):
        ax.text(
            0.82,
            0.96 - j * 0.12,
            rf"{weights[rec]*100:.1f}%",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=12,
            color=colors[j],
        )

plt.xlim(-1,1)
plt.ylim(0,10)
plt.xticks(fontsize=12)
fig.subplots_adjust(hspace=0)
plt.tight_layout()
plt.savefig(plot_file)