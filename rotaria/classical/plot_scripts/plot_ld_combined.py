from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import unphased_local as uph
import sys

# python plot_ld_combined.py silwood
# python plot_ld_combined.py total

mode = sys.argv[1] # silwood or total
pops = ['rm_SIL', 'rs_SIL'] if mode == 'silwood' else ['rm_total', 'rm_SIL', 'rs_total', 'rs_SIL']
targets = [f'{pop}_K_10000_window_100' for pop in pops]
colors = [uph.color_rm if pop.startswith('rm') else uph.color_rs for pop in pops]

def get_data(target):
    data = []
    for f in Path(f"../ld/{target}").glob("*.txt"):
        data.append(np.loadtxt(f))
    data = np.array(data)
    return data

plt.figure(figsize=(6,5))

for target, color in zip(targets, colors):
    data = get_data(target)
    means = np.mean(data, axis=0)
    print(target)
    for i, m in enumerate(means[:20]):
        print(f"{i * 100}-{(i + 1) * 100} bp\t{m:.3f}")
    lower_ci = np.percentile(data, 2.5, axis=0)
    upper_ci = np.percentile(data, 97.5, axis=0)

    plt.plot(range(len(means)), means, linestyle="-" if "SIL" in target else "--", linewidth=2, color="black")
    plt.fill_between(range(len(means)), lower_ci, upper_ci, color=color, alpha=0.3)

ymin, ymax = 0.0, 0.7
plt.ylim(ymin, ymax)
plt.xticks(ticks=range(0, 120, 20), labels=[0, 2, 4, 6, 8, 10], fontsize=18)
plt.yticks(fontsize=18)
# plt.title(f"")
# plt.xlabel("Distance [kb]")
plt.tight_layout()
plt.savefig(f"./combined.{mode}.png")