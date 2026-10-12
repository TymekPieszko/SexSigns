from pathlib import Path
import tskit, msprime, json, sys
import numpy as np
from scipy.optimize import brentq
from multiprocessing import Pool

##############
### About:
# Logic:
# Mean het for phased blocks (min_info=10) in the two species is calculated in the data directory;
# mean value across R. mag. and R. soc. is used as the target.
# Phased blocks with min_info=10 cover approx. 7% and 11% in R. mag. and R. soc., respectively;
# again, 10% is used as an approximate mean between the two.

# This is only run for the GC model,
# the output of which is used for ABC.

##############
### Run as:
# python tune_mut_to_rotaria.py GC 0.011776 5000 0.1 5 16

##############
### Input:
model = sys.argv[1]
target_het = float(sys.argv[2])
window_size = int(sys.argv[3])
bias = float(sys.argv[4])
n_reps = int(sys.argv[5])
n_proc = int(sys.argv[6])
ts_dir = Path(f"sim_output/{model}/1.0.sub_ts/INDS~2")
windows = np.arange(0, 1e6+window_size, window_size)

##############
### Output:
out_file = f"tuned_mut_{model}_target~{target_het}_bias~{bias}.json"

##############
def calc_het(ts_dir, mut):
    reps = []
    for ts_path in Path(ts_dir).glob("*.trees"):
        if int(ts_path.stem) >= n_reps:
            continue
        ts = tskit.load(ts_path)
        ts = msprime.sim_mutations(ts, rate=mut, keep=False)
        hets = ts.diversity([[0,1],[2,3]], windows=windows, mode="site", span_normalise=True)
        hets = np.mean(hets, axis=1)
        top_n = int(len(hets) * bias)
        rep = np.mean(np.sort(hets)[-top_n:])
        reps.append(rep)
    return np.nanmean(reps)

def reach(ts_dir):
    def objective(mut):
        return calc_het(ts_dir, mut) - target_het
    return objective

def worker(ts_dir):
    sex = float(ts_dir.parts[-3].split("~")[1])
    rec = float(ts_dir.parts[-2].split("~")[1])
    target_mut = brentq(reach(ts_dir), 5e-8, 1e-4, xtol=1e-8)
    print(f"SEX: {sex}, REC: {rec}, MUT: {target_mut}", flush=True)
    return f"{sex}_{rec}", target_mut

##############
### Main:
if __name__ == "__main__":
    with Pool(processes=n_proc) as pool:
        total = dict(pool.imap_unordered(worker, sorted(ts_dir.glob("SEX~*/REC~*/MUT~5e-07"))))
    with open(out_file, "w") as f:
        json.dump(total, f, indent=1)
