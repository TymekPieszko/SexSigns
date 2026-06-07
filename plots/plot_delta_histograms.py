import numpy as np
import matplotlib.pyplot as plt
import sys
import json
from pathlib import Path
from sexsigns_utils.plot import sex_labels_in_N, loh_labels_in_N

#############
### ABOUT ###
#############
# Plot the distributions of delta values across the sigma ~ gamma space.

### Run as:
# python plot_delta_histograms.py GC_dbl_b~1.0

###################
### Input:
in_name = sys.argv[1]
in_data = f"../stats/delta/{in_name}.ind_pairs.json"   # adjust if needed
bins = np.linspace(-1,2,31)
print(bins)

###################
### Output:
plot_file = f"./delta/{in_name}.hist.png"
Path(plot_file).parent.mkdir(exist_ok=True)


def plot_delta_grid(in_data, bins):
    keys = list(in_data.keys())
    sex_vals = sorted(set(float(k.split("_")[0]) for k in keys))
    rec_vals = sorted(set(float(k.split("_")[1]) for k in keys))[::-1]

    nrows = len(rec_vals)
    ncols = len(sex_vals)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3 * ncols, 3 * nrows),
        sharex=True,
        sharey=True
    )

    axes = np.atleast_2d(axes)
    for i, rec in enumerate(rec_vals):
        for j, sex in enumerate(sex_vals):
            ax = axes[i, j]
            key = f"{sex}_{rec}"

            y = np.array(in_data[key])
            # print(key)
            # print(np.nanmin(y), np.nanmax(y))
            print(len(y))
            ax.hist(y, bins=bins, color="black", alpha=0.8)
            ax.set_yticks([])
            ax.set_xticks([-1,0,1,2])

            if i == 0:
                ax.set_title(sex_labels_in_N[j], fontsize=32, pad=20)
            if j == 0:
                ax.set_ylabel(loh_labels_in_N[::-1][i], fontsize=32, labelpad=20)
            # if j == len(sex_vals) - 1:
            #     ax.set_ylabel(loh_labels_in_N[::-1][i], fontsize=32, labelpad=20)
            #     ax.yaxis.set_label_position("right")

    fig.subplots_adjust(
        left=0.05,
        right=0.95,
        bottom=0.05,
        top=0.95,
        wspace=0.12,
        hspace=0.12,
    )
    return fig, axes


in_data = json.load(open(in_data))
fig, axes = plot_delta_grid(in_data, bins)

for ax in axes.flat:
    ax.tick_params(axis='both', labelsize=24)

plt.savefig(plot_file, dpi=300, bbox_inches="tight")