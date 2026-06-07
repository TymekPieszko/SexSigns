from pathlib import Path
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
import tskit, msprime, pyslim, sys, json, warnings
from sexsigns_utils.calc import calc_fis, subsample_ts, get_biallelic_01, to_allel

warnings.filterwarnings("ignore")

#############
### ABOUT ###
#############
# Calculate F_IS for individual samples of different sizes (n_lst).

### Run as:
# python calc_fis_sampling.py GC 24
# python calc_fis_sampling.py CO 24

###################
### Input:
model = sys.argv[1]
n_workers = int(sys.argv[2])
target = Path(f"../sim_pipeline_2/sim_output/{model}/0.ts/")
n_lst = [5, 10, 15, 50, 100]

###################
### Output:
out_file = f"./fis_sampling/{model}.json"
Path(out_file).parent.mkdir(exist_ok=True, parents=True)

###################
def worker(args):
    ts_file, n = args

    ts = tskit.load(str(ts_file))
    ts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-6)
    ts = subsample_ts(ts, n)
    ts = msprime.sim_mutations(ts, rate=5e-7)

    genos = ts.genotype_matrix()
    genos = get_biallelic_01(genos, ts.sites_position)[0]
    fis = calc_fis(to_allel(genos))

    return float(fis)

###################
### Main:

total_reps = {}

for n in n_lst:
    for ts_dir in tqdm(list(target.glob("SEX~0.0/REC~*"))):
        rec = float(ts_dir.parts[-1].split("~")[1])
        key = f"{n}_{rec}"

        ts_files = sorted(ts_dir.glob("*.trees"), key=lambda f: float(f.stem))
        tasks = [(f, n) for f in ts_files]

        with Pool(processes=n_workers) as pool:
            reps = list(pool.imap_unordered(worker, tasks))

        total_reps[key] = reps

with open(out_file, "w") as f:
    json.dump(total_reps, f)