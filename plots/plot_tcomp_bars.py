from pathlib import Path
import numpy as np
import json
import argparse
from sexsigns_utils.calc import calc_category_props
from sexsigns_utils.plot import sex_labels_in_N_sparse, loh_labels_in_N_sparse
import matplotlib.pyplot as plt

#############
### ABOUT ###
#############
# Plot tree composition across the sigma ~ gamma space.
# Takes a dictionary with "sex_rec" keys and corresponding lists of replicates as input.

### Run as:
# python plot_tcomp_bars.py --target GC_sim_b~0.1 --diff_from GC_sim_b~1.0
# python plot_tcomp_bars.py --target GC_singer_b~1.0 --diff_from GC_sim_b~1.0
# --diff_from is optional

parser = argparse.ArgumentParser(description="Plot tree composition bars.")
parser.add_argument("--target", required=True)
parser.add_argument("--diff_from", required=False)
args = parser.parse_args()

in_name_1 = args.target
in_name_2 = args.diff_from

###################
### Input:
in_data_1 = (f"../stats/tcomp/{in_name_1}.json")
in_data_1 = json.load(open(in_data_1))

if in_name_2 is not None:
    in_data_2 = (f"../stats/tcomp/{in_name_2}.json")
    in_data_2 = json.load(open(in_data_2))

###################
### Output:
plot_file = Path(f"./tcomp/bars/{in_name_1}.png")
plot_file.parent.mkdir(exist_ok=True, parents=True)

###################
# Get sorted sex_lst & rec_lst;
# necessary for plotting.

sex_lst = sorted({float(k.split("_")[0]) for k in in_data_1.keys()})
rec_lst = sorted({float(k.split("_")[1]) for k in in_data_1.keys()})[::-1]

# Prepare sparse labels
sex_to_label = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
rec_to_label = [1e-6, 1e-7, 1e-8, 1e-9, 0.0]

sex_label_dict = {
    sex: label
    for sex, label in zip(sex_to_label, sex_labels_in_N_sparse)
}

rec_label_dict = {
    rec: label
    for rec, label in zip(rec_to_label, loh_labels_in_N_sparse[::-1])
}

###################
### Plot

colors = {
    "CL": "#377eb8",
    "AR": "#4daf4a",
    "SX": "#ff7f00"
}

fig, axes = plt.subplots(
    nrows=len(rec_lst),
    ncols=len(sex_lst),
    sharex=True,
    sharey=True,
    figsize=(48, 40)
)

for row_idx, rec in enumerate(rec_lst):

    for col_idx, sex in enumerate(sex_lst):

        key = f"{sex}_{rec}"

        ax = axes[row_idx, col_idx]

        # Calc CL-AR-SX proportions (tree composition)
        clarsx_1 = calc_category_props(in_data_1[key])

        # Plot
        bars = ax.bar(
            ["CL", "AR", "SX"],
            clarsx_1.values(),
            color=colors.values()
        )

        ax.set_ylim(0, 1)

        ax.set_yticks([])
        ax.set_xticks([])

        ###################
        ### OPTIONAL:
        ### annotate with difference from ground truth

        if in_name_2 is not None:

            clarsx_2 = calc_category_props(in_data_2[key])

            for bar, v1, v2 in zip(
                bars,
                clarsx_1.values(),
                clarsx_2.values()
            ):

                d = v1 - v2

                height = bar.get_height()

                if height < 0.85:
                    y_pos = height
                    color = "black"
                else:
                    y_pos = height - 0.15
                    color = "white"

                if abs(d) > 0.05:

                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        y_pos,
                        f"{d:.2f}",
                        ha="center",
                        va="bottom",
                        fontsize=40,
                        fontweight="bold",
                        color=color
                    )

        ###################
        ### Add x & y labels!

        if row_idx == len(rec_lst) - 1:

            for k, label in sex_label_dict.items():

                if np.isclose(sex, k):

                    ax.set_xlabel(
                        label,
                        fontsize=60,
                        rotation=0,
                        labelpad=40
                    )

                    break

        if col_idx == 0:

            for k, label in rec_label_dict.items():

                if np.isclose(rec, k):

                    ax.set_ylabel(
                        label,
                        fontsize=60,
                        rotation=0,
                        labelpad=180
                    )

                    break

###################

fig.subplots_adjust(
    left=0.10,   # Space for y labels
    right=0.98,
    top=0.95,
    bottom=0.09,
    wspace=0.0,  # No column gaps
    hspace=0.0,  # No row gaps
)

# Adjust the width of subplot edges
for ax in axes.flat:

    for spine in ax.spines.values():

        spine.set_linewidth(3)

plt.suptitle(
    f"{in_name_1} vs {in_name_2}",
    fontsize=62
)

plt.savefig(plot_file, dpi=240)

plt.close()