from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# python plot_ld.py

target_dir = Path("../LD")

def get_data(target):
    data = []
    for f in Path(target).glob("*.txt"):
        data.append(np.loadtxt(f))
    data = np.array(data)
    return data

for target in target_dir.glob("*"):
    # print(target)
    if not target.is_dir():
        continue
    data = get_data(target)
    means = np.mean(data, axis=0)
    lower_ci = np.percentile(data, 2.5, axis=0)
    upper_ci = np.percentile(data, 97.5, axis=0)
    means, lower_ci, upper_ci = means[:200], lower_ci[:200], upper_ci[:200]
    print(means)
    print(lower_ci)
    print(upper_ci)
    
    # Plot
    plt.figure(figsize=(10,8))

    plt.plot(range(len(means)), means, linewidth=2, color="black")
    plt.fill_between(range(len(means)), lower_ci, upper_ci, color="gray", alpha=0.3)

    first_kb_mean = np.nanmean(means[:10])
    last_kb_mean = np.nanmean(means[-10:])

    plt.text(
        0.98, 0.98,
        f"0-1 kb: {first_kb_mean:.3f}\n"
        f"19-20 kb: {last_kb_mean:.3f}",
        transform=plt.gca().transAxes,
        ha="right",
        va="top",
        fontsize=18,
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none")
    )

    # ymin, ymax = 0.1, 0.7
    plt.ylim(0,1)

    plt.xticks(ticks=range(0,210,10), labels=range(0,21,1), fontsize=26)
    plt.yticks(fontsize=26)

    plt.title(f"{Path(target).name}")

    # plt.xlabel("Distance [kb]")

    plt.tight_layout()
    plt.savefig(target_dir / f"{target.name}.png")