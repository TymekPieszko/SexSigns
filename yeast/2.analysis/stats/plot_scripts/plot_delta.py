from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json

data_dir = Path("../stats")
plot_dir = Path("../delta_plots")

life_cycle_path = Path("../../../config/life_cycle_dict.json")
with open(life_cycle_path, "r") as f:
    life_cycle_lookup = json.load(f)

color_dict = {
    "conventional": "#7F9E94",
    "preferred sexual": "#D26A43",
    "preferred asexual": "#2B2B2B",
    "not measured": "azure",
}

plot_dir.mkdir(parents=True, exist_ok=True)

targets = {
    "ANN", "CQM", "CQE", "CQK", "CQF", "CQN",  # West African Cocoa
    "CNB", "BVG", "CNE", "BVF", "CNS", "BVH"   # Brazilian Bioethanol
}

for file in data_dir.glob("*pair*"):

    try:
        df = pd.read_csv(file, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        continue

    if len(df) < 4:
        continue

    df = df.dropna(subset=["delta"])

    if len(df) == 0:
        continue


    if "Sake" in file.stem:
        bins = np.linspace(-0.75, 1.5, 24)
    else:
        bins = np.linspace(-0.5, 1.5, 21)

    bin_width = bins[1] - bins[0]

    all_inds = pd.unique(df[["ind_A", "ind_B"]].values.ravel())

    lc_counts = pd.Series(
        [life_cycle_lookup.get(ind, "not measured") for ind in all_inds]
    ).value_counts().to_dict()

    def pair_priority(row):

        s1 = row["ind_A"]
        s2 = row["ind_B"]

        special = (s1 in targets) and (s2 in targets)

        lc1 = life_cycle_lookup.get(s1, "not measured")
        lc2 = life_cycle_lookup.get(s2, "not measured")

        c1 = lc_counts.get(lc1, 0)
        c2 = lc_counts.get(lc2, 0)

        rarity = min(c1, c2)

        return (
            0 if special else 1,
            rarity
        )

    priorities = df.apply(pair_priority, axis=1, result_type="expand")
    priorities.columns = ["priority_special", "priority_rarity"]

    df = pd.concat([df, priorities], axis=1)

    df = df.sort_values(
        ["priority_special", "priority_rarity"],
        ascending=[True, True],
        kind="stable"
    ).drop(columns=["priority_special", "priority_rarity"])

    fig, ax = plt.subplots(figsize=(10, 6))

    # track current stack height per bin
    heights = np.zeros(len(bins) - 1, dtype=int)

    # draw contours last
    contours = []

    for _, row in df.iterrows():

        delta = row["delta"]

        s1 = row["ind_A"]
        s2 = row["ind_B"]

        lc1 = life_cycle_lookup.get(s1, "not measured")
        lc2 = life_cycle_lookup.get(s2, "not measured")

        special = (s1 in targets) and (s2 in targets)

        bin_idx = np.digitize(delta, bins) - 1

        if bin_idx < 0 or bin_idx >= len(heights):
            continue

        x_left = bins[bin_idx]
        y_bottom = heights[bin_idx]

        # same life cycle -> solid square
        if lc1 == lc2:

            rect = patches.Rectangle(
                (x_left, y_bottom),
                bin_width,
                1,
                facecolor=color_dict[lc1],
                edgecolor="white",
                linewidth=1.5,
            )

            ax.add_patch(rect)

        # different life cycles -> half-half square
        else:

            count1 = lc_counts.get(lc1, 0)
            count2 = lc_counts.get(lc2, 0)

            colors = [color_dict[lc1], color_dict[lc2]]

            # rarer color on the left
            if count2 < count1:
                colors = colors[::-1]

            # left half
            rect1 = patches.Rectangle(
                (x_left, y_bottom),
                bin_width / 2,
                1,
                facecolor=colors[0],
                edgecolor="none",
                linewidth=0,
            )

            # right half
            rect2 = patches.Rectangle(
                (x_left + bin_width / 2, y_bottom),
                bin_width / 2,
                1,
                facecolor=colors[1],
                edgecolor="none",
                linewidth=0,
            )

            # outer white border only
            border = patches.Rectangle(
                (x_left, y_bottom),
                bin_width,
                1,
                fill=False,
                edgecolor="white",
                linewidth=1.5,
            )

            ax.add_patch(rect1)
            ax.add_patch(rect2)
            ax.add_patch(border)

        # fluorescent contour
        if special:

            contour = patches.Rectangle(
                (x_left, y_bottom),
                bin_width,
                1,
                fill=False,
                edgecolor="#CCFF00",
                linewidth=3,
                zorder=100
            )

            contours.append(contour)

        heights[bin_idx] += 1

    # draw contours on top
    for contour in contours:
        ax.add_patch(contour)

    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, heights.max() + 1)

    # ax.set_xlabel(r"$\it{\Delta}$")
    # ax.set_ylabel("Count")

    ax.set_title(file.stem)

    out_file = plot_dir / f"{file.stem}.png"

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.savefig(out_file, bbox_inches="tight")
    plt.close()