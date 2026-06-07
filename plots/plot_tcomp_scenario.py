from pathlib import Path
import json
import sys
from sexsigns_utils.calc import calc_category_props
import matplotlib.pyplot as plt

#############
### ABOUT ###
#############
# Plot tree composition for a given scenario.

### Run as:
# python plot_tcomp_scenario.py GC_sim_b~1.0 0.3162_1e-06
# python plot_tcomp_scenario.py GC_sim_b~1.0 0.3162_0.0

in_name = sys.argv[1]
key = sys.argv[2]

###################
### Input:
with open(f"../stats/tcomp/{in_name}.json", "r") as f:
    data = json.load(f)

###################
### Output:
plot_file = Path(f"./tcomp/scenario/{key}_{in_name}.png")
plot_file.parent.mkdir(exist_ok=True, parents=True)

###################
### Plot
colors = ["#377eb8", "#4daf4a", "#ff7f00"]
rec_labels = {1e-6:"5.0e-03", 3.162e-07:"1.6e-03", 1e-7:"5.0e-04", 3.162e-08:"1.6e-04", 1e-8:"5.0e-05", 3.162e-09:"1.6e-05", 1e-9:"5.0e-06", 0.0:"0.0"}
fig, ax = plt.subplots(figsize=(9, 10))

clarsx= calc_category_props(data[key])
plt.bar(["CL", "AR", "SX"], clarsx.values(), width=0.8, color=colors)

# X axis
plt.xticks(range(3), ["CL", "AR", "SX"], fontsize=72)
ax.tick_params(axis='x', pad=15)

# Y axis
plt.ylim(0, 1)
plt.yticks(fontsize=60)

plt.title(f"{key}; {in_name}")
plt.tight_layout()
plt.savefig(plot_file)
plt.close()