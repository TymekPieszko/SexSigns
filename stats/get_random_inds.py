import numpy as np
from pathlib import Path
import sys

#############
### ABOUT ###
#############
# Draw random individuals for calculating statistics based on the simulation pipeline.

### Run as:
# python get_random_inds.py 100 30 1000

sim_reps = int(sys.argv[1])
n_ind_pairs = int(sys.argv[2])
total_inds = int(sys.argv[3]) 
out_file = Path(f"inds_{n_ind_pairs*2}.tsv")
rng = np.random.default_rng(42)

if 2 * n_ind_pairs > total_inds:
    raise ValueError("Need at least 2*n_ind_pairs individuals to sample without replacement.")

with open(out_file, "w") as f:
    for _ in range(sim_reps):
        chosen = rng.choice(total_inds, size=2 * n_ind_pairs, replace=False)
        pairs = chosen.reshape(n_ind_pairs, 2)
        pairs.sort(axis=1)  # ensure a <= b within each pair
        f.write("\t".join(map(str, pairs.reshape(-1))) + "\n")