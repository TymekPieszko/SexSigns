from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import sys, json
from rotspec_utils.calc import calc_category_props, rm_silwood, rs_silwood

##############
### Run as:
# python plot_tcomp_rot.py Rmag silwood 10
# python plot_tcomp_rot.py Rsoc silwood 10
# python plot_tcomp_rot.py Rmag total 10
# python plot_tcomp_rot.py Rsoc total 10

##############
### Input:
sp = sys.argv[1]
mode = sys.argv[2] # silwood or total
min_info = int(sys.argv[3])
in_dir = Path(f"../tcomp/{sp}/i~{min_info}")

##############
### Output:
plot_1 = Path(f"../plots/tcomp/{sp}.{mode}.i~{min_info}.TOTAL.png")
plot_2 = Path(f"../plots/tcomp/{sp}.{mode}.i~{min_info}.PAIR.png")
plot_1.parent.mkdir(exist_ok=True)

##############
### Compute props_lst and props_total
tcomp_lst = []
pair_lst = []
for file in sorted(Path(in_dir).glob("*.json")):
    with open(file, "r") as f:
        tcomp = json.load(f)
    tcomp_lst.append(tcomp)
    pair_lst.append(file.stem)
silwood = {"Rmag": rm_silwood, "Rsoc": rs_silwood}[sp]
is_silwood = np.array([set(pair.split("~")) <= set(silwood) for pair in pair_lst])
if mode == "silwood":
    tcomp_lst = [t for t, s in zip(tcomp_lst, is_silwood) if s]
    pair_lst = [p for p, s in zip(pair_lst, is_silwood) if s]
    is_silwood = is_silwood[is_silwood]
props_lst = [calc_category_props(tcomp) for tcomp in tcomp_lst]
props_total = {k: sum(d[k] for d in props_lst) for k in props_lst[0]}

##############
### Categories, as a (n_pairs, 3) array
cats = ["CL", "AR", "SX"]
colors = ["#377eb8", "#4daf4a", "#ff7f00"]
props_arr = np.array([[p[c] for c in cats] for p in props_lst], dtype=float)

##############
### Plot 1: total
# One band per category: min-max across pairs, with the mean dashed.
mins, maxs, means = props_arr.min(axis=0), props_arr.max(axis=0), props_arr.mean(axis=0)
fig = plt.figure(figsize=(6.4, 6.9))
ax = fig.add_axes([1.5 / 6.4, 0.4 / 6.9, 4.5 / 6.4, 6 / 6.9])
for i in range(len(cats)):
    x = [i, i + 1]
    ax.fill_between(x=x, y1=mins[i], y2=maxs[i], color=colors[i], alpha=0.2, linewidth=0)
    ax.plot(x, [mins[i], mins[i]], color=colors[i], linewidth=6)
    ax.plot(x, [maxs[i], maxs[i]], color=colors[i], linewidth=6)
    ax.plot(x, [means[i], means[i]], color=colors[i], linestyle=(0, (1, 0.5)), linewidth=14)
ax.set_xticks(np.arange(len(cats)) + 0.5)
ax.tick_params(axis="x", labelbottom=False, length=18, width=3)
ax.tick_params(axis="y", labelsize=38)
ax.set_xlim(0, 3); ax.set_ylim(0, 1); ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_title(f"{sp} i~{min_info} (n={len(props_arr)})", pad=15)
fig.savefig(plot_1, dpi=150)

### Save means alongside the plot
with open(plot_1.with_suffix(".means.json"), "w") as f:
    json.dump(dict(zip(cats, means.tolist())), f, indent=1)

##############
### Plot 2: pairs
# One stacked horizontal bar per pair; the three props sum to 1.
y = np.arange(len(pair_lst))
fig, ax = plt.subplots(figsize=(10, 0.32 * len(pair_lst) + 2))
left = np.zeros(len(pair_lst))
for i, cat in enumerate(cats):
    ax.barh(y, props_arr[:, i], left=left, color=colors[i], label=cat, height=0.8)
    left += props_arr[:, i]
ax.set_yticks(y); ax.set_yticklabels(pair_lst, fontsize=10)
ax.invert_yaxis()          # first pair at the top
ax.set_xlim(0, 1)
ax.set_xlabel("proportion", fontsize=12)
ax.legend(ncol=len(cats), fontsize=12, loc="lower center",
          bbox_to_anchor=(0.5, 1.01), frameon=False)
ax.set_title(f"{sp} i~{min_info} (n={len(props_arr)})", pad=40)
fig.tight_layout()
fig.savefig(plot_2, dpi=150)

print(f"{plot_1}\n{plot_2}", flush=True)