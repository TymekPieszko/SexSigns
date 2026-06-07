from pathlib import Path
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
import tskit, msprime, pyslim, json, sys, warnings
from sexsigns_utils.calc import read_inds, subsample_to_inds, get_biallelic_01, calc_delta

warnings.simplefilter('ignore')

#############
### ABOUT ###
#############
# Calculate delta for target individuals (grouped into pairs) drawn from simulation replicates.
# 'interval_type' specifies how intervals are defined (e.g., ARG- or SNP-based). 

### Run as:
# python calc_delta.py GC dbl 1.0 16
# python calc_delta.py GC dbl_trp 1.0 16
# python calc_delta.py GC CL 1.0 16
# python calc_delta.py GC CL_SX 1.0 16

###################
### Input:
sim_model = sys.argv[1]
interval_type = sys.argv[2]
bias = float(sys.argv[3])
n_workers = int(sys.argv[4])
inds_file = "./inds_60.tsv"
target = Path(f"../sim_pipeline_2/sim_output/{sim_model}/0.ts/")
ordered_window_dict = f"./ordered_windows/ordered_windows_{sim_model}.ts.json"
with open(ordered_window_dict, "r") as f:
    ordered_window_dict = json.load(f)
delta_interval_dict = f"./delta_intervals/{interval_type}_intervals_{sim_model}.json"
with open(delta_interval_dict, "r") as f:
    delta_interval_dict = json.load(f)

###################
### Output:
out_file_1 = Path(f"./delta/{sim_model}_{interval_type}_b~{bias}.json") # Mean per rep
out_file_2 = Path(f"./delta/{sim_model}_{interval_type}_b~{bias}.ind_pairs.json") # All values per rep
out_file_1.parent.mkdir(exist_ok=True)

###################
def worker(args):
    ts, key = args
    i = ts.stem
    ordered_windows = ordered_window_dict[key][i]
    delta_intervals = delta_interval_dict[key][i]

    ts = tskit.load(ts)
    ts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-6)
    ts = msprime.sim_mutations(ts, rate=5e-7)

    delta_j_lst = []
    for j, (ind_1, ind_2) in enumerate(ind_pairs[int(i)]):
        ts_2 = subsample_to_inds(ts, ind_1, ind_2)

        ow_j = ordered_windows[str(j)]
        di_j = delta_intervals[str(j)]

        n_keep = max(1, round(len(ow_j) * bias))
        windows = sorted(ow_j[:n_keep], key=lambda x: x[0])  # coordinate sorting for keep_intervals
        ts_2 = ts_2.keep_intervals(windows).trim()

        genos = ts_2.genotype_matrix()
        pos = ts_2.sites_position
        genos, pos = get_biallelic_01(genos, pos)

        if len(genos) == 0:
            delta_j_lst.append(np.nan)
            continue

        delta_j = calc_delta(genos, pos, di_j)[0]
        delta_j_lst.append(delta_j)

    return i, np.nanmean(delta_j_lst), delta_j_lst

###################
### Main:

ind_pairs = read_inds(inds_file)
means_total = {}
pairs_total = {}

for dir in target.glob("SEX~*/REC~*"):
    sex = float(dir.parts[-2].split("~")[1])
    rec = float(dir.parts[-1].split("~")[1])
    key = f"{sex}_{rec}"
    print(key)

    ts_files = sorted(dir.glob("*.trees"), key=lambda f: int(f.stem))

    means_lst = []
    pairs_lst = []

    tasks = [(f, key) for f in ts_files]
    with Pool(processes=n_workers) as pool:
        for i, mean_delta, pair_deltas in tqdm(
            pool.imap_unordered(worker, tasks),
            total=len(tasks),
        ):
            means_lst.append(mean_delta)
            pairs_lst.extend(pair_deltas)

    means_total[key] = means_lst
    pairs_total[key] = pairs_lst

with open(out_file_1, "w") as f:
    json.dump(means_total, f)

with open(out_file_2, "w") as f:
    json.dump(pairs_total, f)