from pathlib import Path
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
import tskit, msprime, pyslim, sys, json, warnings
from sexsigns_utils.calc import get_biallelic_01, subsample_ts, to_allel, calc_r2_windowed

warnings.filterwarnings("ignore")

#############
### ABOUT ###
#############
# Calculate mean LD decay curves for simulation scenarios.

### Run as:
# python calc_ld_curves.py GC 16

###################
### Input:
model = sys.argv[1]
n_workers = int(sys.argv[2])
inds_file = "./inds_60.tsv"
target = Path(f"../sim_pipeline_2/sim_output/{model}/0.ts/")
n = 10 # Number of inds
k = 100 # Number of SNPs
window_len = 1000
window_num = 1001

###################
### Output:
out_file = f"./ld_curves/{model}.json"
Path(out_file).parent.mkdir(exist_ok=True, parents=True)

###################
def worker(ts_file):

    ts = tskit.load(str(ts_file))
    ts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-6)
    ts = subsample_ts(ts, n)
    ts = msprime.sim_mutations(ts, rate=5e-7)
    genos = ts.genotype_matrix()
    pos = ts.sites_position
    genos, pos = get_biallelic_01(genos, pos)

    r2 = calc_r2_windowed(to_allel(genos), pos, k, window_len, window_num)
    return r2
###################

###################
### Main:

total_reps = {}
for ts_dir in target.glob("SEX~*/REC~*"):
    sex = float(ts_dir.parts[-2].split("~")[1])
    rec = float(ts_dir.parts[-1].split("~")[1])
    key = f"{sex}_{rec}"
    print(key)

    ts_files = sorted(ts_dir.glob("*.trees"), key=lambda f: float(f.stem))

    reps = []
    with Pool(processes=n_workers) as pool:
        for val in tqdm(
            pool.imap_unordered(worker, ts_files),
            total=len(ts_files),
        ):
            reps.append(val)

    total_reps[key] = np.nanmean(reps, axis=0).tolist()
    # print(np.nanmean(reps, axis=0))

with open(out_file, "w") as f:
    json.dump(total_reps, f)
