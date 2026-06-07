from pathlib import Path
import sys, json
import numpy as np
import matplotlib.pyplot as plt
from functions_val import *
from sexsigns_utils.plot import loh_labels_in_N

# python plot_validation.py GC Hi_windows_GC_mut_5e-07_inds_100
# python plot_validation.py CO Hi_windows_CO_mut_5e-07_inds_100
model = sys.argv[1]
in_name = sys.argv[2]
L = 1e6
TRACT = 5000
MUT = 5e-07
window_length = 10000

sim_data = f"./out/{in_name}.txt"
plot_file = Path(f"./out/{in_name}.png")

with open(sim_data, "r") as f:
    sim_data = json.load(f)

sim_data = {
    SEX: {REC: np.nanmean(reps, axis=0) for REC, reps in REC_dict.items()}
    for SEX, REC_dict in sim_data.items()
}

SEX_lst = ["0.0"]
REC_lst = list(next(iter(sim_data.values())).keys())
combos = [(s, r) for s in SEX_lst for r in REC_lst]

fig, ax = plt.subplots(
    nrows=2, ncols=int(len(REC_lst) / 2), sharex=False, figsize=(15, 9)
)
ax = ax.flatten()

for i, combo in enumerate(combos):
    print(combo)
    gamma_label = loh_labels_in_N[REC_lst.index(combo[1])]

    if model == "GC":
        positions = np.arange(0, L, window_length, dtype=int)
        hi_analytical = calc_hi(MUT, calc_gamma_GC(float(combo[1]), TRACT))
        hi_simulated = sim_data[str(combo[0])][str(combo[1])]
        ax[i].plot(
            positions,
            hi_simulated,
            color="cornflowerblue",
            marker="o",
            markersize=3,
            linestyle="",
        )
        ax[i].axhline(y=hi_analytical, color="black", linestyle="-")
        ax[i].set_title(r"$\gamma$ = " + str(gamma_label), fontsize=14)
        ax[i].set_ylabel(r"$H_{I}$", fontsize=12)
        ax[i].set_xlabel("Position (bp)", fontsize=12, labelpad=0.5)
        ax[i].set_ylim(bottom=0)

    elif model == "CO":
        positions = np.arange(
            window_length, L + window_length, window_length, dtype=int
        ) - (window_length / 2)
        print(len(positions))
        hi_analytical = [
            calc_hi(MUT, calc_gamma_CO_NCI(float(combo[1]), pos)) for pos in positions
        ]
        print(len(hi_analytical))
        hi_simulated = sim_data[str(combo[0])][str(combo[1])]
        print(len(hi_simulated))
        ax[i].plot(
            positions,
            hi_simulated,
            color="cornflowerblue",
            marker="o",
            markersize=3,
            linestyle="",
        )
        ax[i].plot(positions, hi_analytical, color="black", linestyle="-")
        ax[i].set_xscale("log")
        ax[i].set_xlim(left=1e4, right=1e6)
        ax[i].set_ylim(bottom=0)
        ax[i].set_title(r"$\gamma$ = " + str(gamma_label), fontsize=14)
        ax[i].set_ylabel(r"$H_{I}$", fontsize=12)
        ax[i].set_xlabel("Position (bp)", fontsize=12, labelpad=0.5)

plt.tight_layout()
plt.savefig(plot_file, dpi=320)
plt.close()