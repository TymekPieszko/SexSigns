from pathlib import Path
import numpy as np
from tqdm import tqdm
from itertools import groupby
from multiprocessing import Pool
import tskit, msprime, pyslim, json, sys, warnings
from sexsigns_utils.calc import read_inds, subsample_to_inds, get_biallelic_01, get_tree_intervals, get_snp_intervals

warnings.simplefilter('ignore')

#############
### ABOUT ###
#############
# Find intervals for computing delta;
# used in calc_delta.py.

### Run as:
# python get_delta_intervals.py GC CL 16
# python get_delta_intervals.py GC CL_SX 16
# python get_delta_intervals.py GC dbl 16
# python get_delta_intervals.py GC dbl_trp 16

###################
### Input:
sim_model = sys.argv[1]
interval_type = sys.argv[2]
n_workers = int(sys.argv[3])
inds_file = "./inds.tsv"
target_dir = Path(f"../sim_pipeline_2/sim_output/{sim_model}/0.ts/")

###################
### Output:
out_file = f"./delta_intervals/{interval_type}_intervals_{sim_model}.json"
Path(out_file).parent.mkdir(exist_ok=True)

###################
def worker(args):
    ts_path, inds = args
    i = int(ts_path.stem)

    ts = tskit.load(ts_path)
    ts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-6)
    ts = msprime.sim_mutations(ts, rate=5e-7)

    out = {}
    for j, (ind_1, ind_2) in enumerate(inds):
        ts_2 = subsample_to_inds(ts, ind_1, ind_2)
        if interval_type == "CL":
            intervals = get_tree_intervals(ts_2, include_SX=False)
        elif interval_type == "CL_SX":
            intervals = get_tree_intervals(ts_2, include_SX=True)
        elif interval_type == "dbl":
            genos = ts_2.genotype_matrix()
            pos = ts_2.sites_position
            seq_len = ts_2.sequence_length
            genos, pos = get_biallelic_01(genos, pos)
            intervals = get_snp_intervals(genos, pos, seq_len, include_trp=False)
        elif interval_type == "dbl_trp":
            genos = ts_2.genotype_matrix()
            pos = ts_2.sites_position
            seq_len = ts_2.sequence_length
            genos, pos = get_biallelic_01(genos, pos)
            intervals = get_snp_intervals(genos, pos, seq_len, include_trp=True)
        out[j] = intervals

    return i, out
###################

###################
### Main:

ind_pairs = read_inds(inds_file)

total = {}
for dir in target_dir.glob("SEX~*/REC~*"):
    sex = float(dir.parts[-2].split("~")[1])
    rec = float(dir.parts[-1].split("~")[1])
    key = f"{sex}_{rec}"
    print(key)
    
    total[key] = {}
    ts_lst = sorted(dir.glob("*.trees"), key=lambda ts: int(ts.stem))
    tasks = [(p, ind_pairs[int(p.stem)]) for p in ts_lst]

    with Pool(processes=n_workers) as pool:
        for rep_id, out in tqdm(
            pool.imap_unordered(worker, tasks),
            total=len(tasks),
        ):
            total[key][rep_id] = out

with open(out_file, "w") as f:
    json.dump(total, f)