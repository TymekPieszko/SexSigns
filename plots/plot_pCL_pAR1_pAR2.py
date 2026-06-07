from pathlib import Path
import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from sexsigns_utils.plot import loh_labels_in_N

#############
### ABOUT ###
#############
# Plot p(CL), p(AR1) and p(AR2) under obligate asexuality.

### Run as:
# python plot_pCL_pAR1_pAR2.py GC_sim_b~1.0

def calc_pCL_pAR1_pAR2(tcomp):
    CL = sum(
        v for k, v in tcomp.items()
        if k in {"b0213", "b1302", "b0312", "b1203", "b5_eq", "b6_eq"}
    )

    AR1 = sum(
        v for k, v in tcomp.items()
        if k in {"u0123", "u0132", "u2301", "u2310"}
    )

    AR2 = sum(
        v for k, v in tcomp.items()
        if k in {"b0123", "b2301", "b4_eq"}
    )

    total = CL + AR1 + AR2
    return CL / total, AR1 / total, AR2 / total

in_name = sys.argv[1]

###################
### Input
with open(f"../stats/tcomp/{in_name}.json") as f:
    data = json.load(f)

###################
### Output
plot_file = Path(f"./tcomp/obligate_asex/{in_name}.pCL_pAR1_pAR2.png")
plot_file.parent.mkdir(parents=True, exist_ok=True)

###################
### Get sorted rec_lst
rec_lst = sorted({float(k.split("_")[1]) for k in data})

###################
### Collect p(CL), p(AR1), p(AR2)
p_CL_lst = []
p_AR1_lst = []
p_AR2_lst = []

for rec in rec_lst:
    key = f"0.0_{rec}"

    p_CL, p_AR1, p_AR2 = calc_pCL_pAR1_pAR2(data[key])

    p_CL_lst.append(p_CL)
    p_AR1_lst.append(p_AR1)
    p_AR2_lst.append(p_AR2)
print(np.max(p_AR1_lst))
###################
### Plot
x = range(len(rec_lst))

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(
    x, p_CL_lst,
    "-o",
    lw=3,
    ms=8,
    color="blue",
    label=r"$p(CL)$"
)

ax.plot(
    x, p_AR1_lst,
    "-o",
    lw=3,
    ms=8,
    color="limegreen",
    label=r"$p(AR_1)$"
)

ax.plot(
    x, p_AR2_lst,
    "-o",
    lw=3,
    ms=8,
    color="darkgreen",
    label=r"$p(AR_2)$"
)

ax.set_xticks(list(x))
ax.set_xticklabels(loh_labels_in_N, rotation=45)

ax.set_ylim(-0.05, 1.05)
ax.set_ylabel("Proportion")
ax.legend()

plt.tight_layout()
plt.savefig(plot_file, dpi=300)
plt.close()