from pathlib import Path
import numpy as np

### Run as:
# python make_inds.py

inds_2 = np.loadtxt("sim_output/inds_pipeline.txt", dtype=int, ndmin=2)
np.savetxt("sim_output/inds_2.txt", inds_2, fmt="%d", delimiter="\t")
inds_10 = []
for rep, targets in enumerate(inds_2):
    rng = np.random.default_rng(rep)
    inds_10.append(np.concatenate([targets, rng.choice(np.setdiff1d(np.arange(1000), targets), 8, replace=False)]))
np.savetxt("sim_output/inds_10.txt", np.array(inds_10), fmt="%d", delimiter="\t")
