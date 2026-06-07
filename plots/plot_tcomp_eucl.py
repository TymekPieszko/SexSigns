from pathlib import Path
import json
import sys
import numpy as np
from scipy.spatial.distance import euclidean
from sexsigns_utils.calc import calc_category_props
from sexsigns_utils.plot import plot_heatmap

#############
### ABOUT ###
#############
# Plot Euclidean distance between the truth (1st arg) and a target (2nd arg).
# The distance is calculated for CL-AR-SX proportions of trees.

### Run as:
# python plot_tcomp_eucl.py GC_sim_b~1.0 GC_sticcs_b~1.0

in_name_1 = sys.argv[1]
in_name_2 = sys.argv[2]

###################
### Input:
data_1 = (f"../stats/tcomp/{in_name_1}.json")
with open(data_1, "r") as f:
    data_1 = json.load(f)
data_2 = (f"../stats/tcomp/{in_name_2}.json")
with open(data_2, "r") as f:
    data_2 = json.load(f)
###################
### Output:
plot_file = Path(f"./tcomp/eucl/{in_name_2}_eucl.png")
plot_file.parent.mkdir(exist_ok=True, parents=True)

###################
### Get sorted sex_lst & rec_lst
sex_lst = sorted({float(k.split("_")[0]) for k in data_1.keys()})
rec_lst = sorted({float(k.split("_")[1]) for k in data_1.keys()})[::-1]

###################
### Compute distance matrix D
D = np.zeros((len(rec_lst), len(sex_lst)))
for i, rec in enumerate(rec_lst):
    for j, sex in enumerate(sex_lst):

        key = f"{sex}_{rec}"

        clarsx_1 = list(calc_category_props(data_1[key]).values())
        clarsx_2 = list(calc_category_props(data_2[key]).values())

        D[i, j] = euclidean(clarsx_1, clarsx_2)

###################
### Plot D
fig, ax = plot_heatmap(
    data=D,
    cmap="viridis",
    vmin=0,
    vmax=0.21,
    xlabel=r"$\sigma$",
    ylabel=r"$\gamma$",
    tick_labels_in_N=True,
    title=f"{in_name_1} vs {in_name_2}",
    title_font=18,
    label_font=62,
    tick_font=32,
    cbar=True, # Legend
    # annot=True,
    annot=False,
)

fig.savefig(plot_file)
