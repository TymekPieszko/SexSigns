from pathlib import Path
import json
import sys
from sexsigns_utils.calc import calc_category_props
from sexsigns_utils.plot import sex_labels_in_N_sparse, loh_labels_in_N_sparse
import matplotlib.pyplot as plt

#############
### ABOUT ###
#############
# Plot a vertical slice throught the sigma ~ gamma space (for sigma = x).
# Optionally plot ground truth.

### Run as:
# python plot_tcomp_vertical_slice.py GC_sim_b~0.1 GC_sim_b~1.0 0.0 True
# python plot_tcomp_vertical_slice.py GC_singer_b~1.0 GC_sim_b~1.0 0.0 True
# python plot_tcomp_vertical_slice.py GC_sticcs_b~1.0 GC_sim_b~1.0 0.0 True

###################
### Input:
if len(sys.argv) == 4:
    in_name_1 = sys.argv[1]
    in_name_2 = None
    sex = float(sys.argv[2])
    annot_SX = sys.argv[3]
elif len(sys.argv) == 5:
    in_name_1 = sys.argv[1]
    in_name_2 = sys.argv[2] # Optional!
    sex = float(sys.argv[3])
    annot_SX = sys.argv[4]
in_data_1 = f"../stats/tcomp/{in_name_1}.json"
in_data_1 = json.load(open(in_data_1))
if in_name_2 is not None:
    in_data_2 = f"../stats/tcomp/{in_name_2}.json"
    in_data_2 = json.load(open(in_data_2))
###################
### Output:
plot_file = Path(f"./tcomp/vertical_slice/{in_name_1}_sex={sex}.png")
plot_file.parent.mkdir(exist_ok=True, parents=True)

###################
# Get sorted sex_lst & rec_lst;
# necessary for plotting.
sex_lst = sorted({float(k.split("_")[0]) for k in in_data_1.keys()})
rec_lst = sorted({float(k.split("_")[1]) for k in in_data_1.keys()})[::-1]
# Prepare labels
sex_to_label = [0.0, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
rec_to_label = [1e-6, 1e-7, 1e-8, 1e-9, 0.0]
sex_label_dict = {sex:label for sex,label in zip(sex_to_label,sex_labels_in_N_sparse)}
rec_label_dict = {rec:label for rec,label in zip(rec_to_label,loh_labels_in_N_sparse[::-1])}

###################
### Plot
colors = {"CL": "#377eb8", "AR": "#4daf4a", "SX": "#ff7f00"}
fig, axes = plt.subplots(nrows=len(rec_lst), ncols=1, sharey=True, figsize=(10, 40))
for i, rec in enumerate(rec_lst):
    ax = axes[i]
    key = f"{sex}_{rec}"

    # Plot target data
    clarsx_1 = calc_category_props(in_data_1[key])
    bars = ax.bar(["CL", "AR", "SX"], clarsx_1.values(), width=0.8, color=colors.values())
    # OPTIONAL: annotate the height of SX bar
    if annot_SX:
        bar = bars[-1]
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,
            height,
            f"{height:.2f}",
            ha="center",
            va="bottom",
            fontsize=48,
            fontweight="bold"
        )
    # OPTONAL: plot ground truth
    if in_name_2 is not None:
        clarsx_2 = calc_category_props(in_data_2[key])
        idx = 0
        for tree_cat, prop_1, prop_2 in zip(clarsx_2.keys(), clarsx_1.values(), clarsx_2.values()):
            if prop_1 > prop_2:
                color="white"
            else:
                color=colors[tree_cat]
            ax.hlines(y=prop_2, xmin=idx-0.4, xmax=idx+0.4, color=color, linestyle=(0,(1,1)), linewidth=10)
            idx += 1
    
    # Y axis
    ax.set_ylim(0, 1)
    # ax.set_yticks(ticks=[0,0.2,0.4,0.6,0.8,1.0], labels=[0,0.2,0.4,0.6,0.8,1.0], fontsize=32)
    ax.set_yticks([]); ax.set_xticks([])

    # X axis
    ax.tick_params(axis="x", labelsize=28)
    if rec in rec_label_dict: 
        ax.set_ylabel(rec_label_dict[rec], fontsize=60, rotation=0, labelpad=180)
    if i == 7:
        ax.set_xlabel(sex_label_dict[sex], fontsize=60, rotation=0, labelpad=40)

fig.subplots_adjust(
    left=0.5,   # Space for y labels
    right=0.88,
    top=0.95,
    bottom=0.09,
    wspace=0.0,  # No column gaps
    hspace=0.0,   # No row gaps
)

# Adjust the width of subplot edges
for ax in axes.flat:
    for spine in ax.spines.values():
        spine.set_linewidth(3)

plt.suptitle(in_name_1, y=0.995)
plt.savefig(plot_file)
plt.close()