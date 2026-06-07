from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# python plot_R2.py abc_R2_rec_0e+00_to_1e-06 abc_R2_rec_3.16e-08_to_3.16e-07 abc_R2_rec_1e-07_to_1e-07

###################
# Input
targets = sys.argv[1:]

###################
# Output
plot_file = "./out/R2_plot.png"

###################
# Main

pos = [0, 1, 2, 3.5, 4.5, 5.5, 7, 8, 9]
bar_width = 1.0

fig, ax = plt.subplots(figsize=(12, 7))

alpha_levels = np.linspace(0.15, 0.35, len(targets))
colors = ["aliceblue", "cornflowerblue", "navy"]
# hatches = ["/", "//", "|"]

for i, (alpha, target) in enumerate(zip(alpha_levels, targets)):

    data = f"./out/abc_results/{target}.tsv"
    df = pd.read_csv(data, sep="\t")

    values = df.iloc[0].values
    color = colors[i]

    bars = ax.bar(
        pos,
        values,
        width=bar_width,
        color="black",
        # hatch=hatches[i],
        edgecolor="black",
        linewidth=2,
        alpha=alpha,
    )

    for j, (x, y) in enumerate(zip(pos, values)):
        ax.hlines(
            y=y,
            xmin=x - bar_width / 2,
            xmax=x + bar_width / 2,
            colors=color,
            linewidth=3,
            label=target if j == 0 else None  # only once for legend
        )

ax.set_xticks(pos)
ax.set_xticklabels(df.columns, rotation=66, fontsize=14)
ax.set_ylim(0, 1)
ax.tick_params(axis="y", labelsize=14)
ax.legend(frameon=False)

plt.tight_layout()
plt.savefig(plot_file, dpi=320)
plt.show()