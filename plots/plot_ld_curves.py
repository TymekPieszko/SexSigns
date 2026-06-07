import numpy as np
import matplotlib.pyplot as plt
import sys
import json
from pathlib import Path
from sexsigns_utils.plot import sex_labels_in_N, loh_labels_in_N

#############
### ABOUT ###
#############
# Plot LD decay curves across the sigma ~ gamma space.

### Run as:
# python plot_ld_curves.py GC

###################
### Input:
model = sys.argv[1]
in_data = f"../stats/ld_curves/{model}.json"
window_len = 1000
window_num = 101

###################
### Output:
plot_file = f"./ld_curves/{model}.png"
Path(plot_file).parent.mkdir(exist_ok=True)


def plot_ld_grid(in_data, window_len):
    keys = list(in_data.keys())
    sex_vals = sorted(set(float(k.split("_")[0]) for k in keys))
    rec_vals = sorted(set(float(k.split("_")[1]) for k in keys))[::-1]
    # print(sex_vals)
    nrows = len(rec_vals)
    ncols = len(sex_vals)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(3 * ncols, 4 * nrows),
        sharex=True, sharey=True
    )

    axes = np.atleast_2d(axes)

    for i, rec in enumerate(rec_vals):
        for j, sex in enumerate(sex_vals):
            ax = axes[i, j]
            key = f"{sex}_{rec}"

            if key not in in_data:
                ax.axis("off")
                continue

            y = np.array(in_data[key])
            x = np.arange(len(y)) * window_len + window_len / 2

            ax.plot(x[:window_num], y[:window_num], color="black")
            # Highlight 1st & 20th windows
            # 1st window (index 0)
            ax.axvspan(0, window_len, color="gold")
            ax.axvspan(19 * window_len, 20 * window_len, color="gold")
            # Annotate with fold decay
            decay_1 = y[0] - y[19]
            ax.text(
                0.98, 0.98,
                # rf"20 kb: {decay_1:.3f}",
                rf"{decay_1:.3f}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=24,
                # color="gold"
            )
            # decay_2 = y[0] - y[499]
            # ax.text(
            #     0.98, 0.86,
            #     rf"500 kb: {decay_2:.3f}",
            #     transform=ax.transAxes,
            #     ha="right",
            #     va="top",
            #     fontsize=18,
            #     # color="gold"
            # )

            if i == 0:
                try:
                    ax.set_title(sex_labels_in_N[j], fontsize=32, pad=20)
                except:
                    print(j)
            if j == len(sex_vals)-1:
                ax.set_ylabel(loh_labels_in_N[::-1][i], fontsize=32, labelpad=20)
                ax.yaxis.set_label_position("right")

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
fig, axes = plot_ld_grid(in_data, window_len)
for ax in axes.flat:
    # ax.set_xticks(np.arange(101, step=20)*window_len)
    # ax.set_xticklabels(np.arange(101, step=20))
    ax.set_xticks([0,50000,100000])
    ax.set_xticklabels([0,50,100])
    ax.tick_params(axis='both', labelsize=24)
plt.savefig(plot_file, dpi=300)