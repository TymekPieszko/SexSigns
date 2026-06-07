from pathlib import Path
import json
import sys
import numpy as np
from sexsigns_utils.plot import plot_heatmap

#############
### ABOUT ###
#############
# Plot a standardised heatmap acros the sex ~ rec parameter space.
# Takes a dictionary with "sex_rec" keys and corresponding lists of replicated values as input.

### Run as:
# python plot_heatmap.py fis GC
# python plot_heatmap.py hi GC
# python plot_heatmap.py LD_0-1000 GC
# python plot_heatmap.py delta GC_CL_b~1.0
# python plot_heatmap.py delta GC_CL_b~0.1
# python plot_heatmap.py kappa GC

###################
### Input:
stat = sys.argv[1]
in_name = sys.argv[2]
data = (f"../stats/{stat}/{in_name}.json")
# data = "/data/biol-bdelloids/scro4331/SexSigns_2025/stats/LD_0-1000/LD_0-1000_GC_mut_5e-07_inds_100.json"
with open(data, "r") as f:
    data = json.load(f)

###################
### Output:
plot_file = Path(f"./{stat}/{in_name}.png")
plot_file.parent.mkdir(exist_ok=True, parents=True)

###################
### Get sorted sex_lst & rec_lst
sex_lst = sorted({float(k.split("_")[0]) for k in data.keys()})
sex_lst = [sex for sex in sex_lst if sex != 1.0]
rec_lst = sorted({float(k.split("_")[1]) for k in data.keys()})[::-1]

###################
### Compute mean values per sex~rec
to_plot = np.zeros((len(rec_lst), len(sex_lst)))
for i, rec in enumerate(rec_lst):
    for j, sex in enumerate(sex_lst):
        key = f"{sex}_{rec}"
        # print(key, np.nanmean(data[key]))
        to_plot[i, j] = np.nanmean(data[key])

# print(to_plot[to_plot>1000])
to_plot[to_plot>1000] = np.nan
# print(to_plot)
print(to_plot.min())
print(to_plot.max())
###################
### Plot
fig, ax = plot_heatmap(
    data=to_plot,
    cmap="plasma",
    vmin=np.nanmin(to_plot),
    vmax=np.nanmax(to_plot),
    # vmax=0.4,
    xlabel=r"$\sigma$",
    ylabel=r"$\gamma$",
    tick_labels_in_N=True,
    title=f"{stat}; {in_name}",
    #########
    # For paper
    title_font=18,
    label_font=62,
    tick_font=32,
    #########
    # For talk
    # title_font=26,
    # label_font=42,
    # tick_font=32,
    #########
    cbar=True, # Legend
    # cbar=False, # Legend
    # annot=True,
    annot=False,
)

fig.savefig(plot_file)
