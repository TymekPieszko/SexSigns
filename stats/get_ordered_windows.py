from pathlib import Path
import numpy as np
from tqdm import tqdm
import tskit, msprime, pyslim, json, sys, warnings
from sexsigns_utils.calc import read_inds

warnings.simplefilter('ignore')

#############
### ABOUT ###
#############
# Compute windows used to implement the heterozygosity bias;
# the windows are returned in order of decreasing heterozygosity.

### Run as:
# python get_ordered_windows.py GC 5000 ts
# python get_ordered_windows.py GC 5000 sub_ts

###################
### Input:
sim_model = sys.argv[1]
window_length = int(sys.argv[2])
mode = sys.argv[3]
if mode == "ts":
    target = "0.ts"
    inds_file = "./inds.tsv"
    glob_pattern = "SEX~*/REC~*"
    sex_idx, rec_idx = -2, -1
elif mode == "sub_ts":
    target = "1.0.sub_ts"
    glob_pattern = "SEX~*/REC~*/MUT~*"
    sex_idx, rec_idx = -3, -2
target_dir = Path(f"../sim_pipeline_2/sim_output/{sim_model}/{target}/")

###################
### Output:
out_file = f"./ordered_windows/ordered_windows_{sim_model}.{mode}.json"
Path(out_file).parent.mkdir(exist_ok=True)

###################
def get_ordered_windows(ts, window_length, A, B):
    # A & B are individual IDs
    breaks = np.arange(0, ts.sequence_length + window_length, window_length)
    nodes_A = ts.individual(A).nodes
    nodes_B = ts.individual(B).nodes
    het = ts.diversity([nodes_A, nodes_B], windows=breaks, mode="branch", span_normalise=True)
    het = np.mean(het, axis=1)
    windows = np.column_stack([breaks[:-1], breaks[1:]])
    idx = np.argsort(het)[::-1]  # Order of decreasing het
    return windows[idx]
###################

###################
### Main:

ind_pairs = read_inds(inds_file) if target == "plot" else None
total = {}
for dir in target_dir.glob(glob_pattern):
    sex = float(dir.parts[sex_idx].split("~")[1])
    rec = float(dir.parts[rec_idx].split("~")[1])
    key = f"{sex}_{rec}"
    print(key)

    total[key] = {}
    ts_lst = sorted(dir.glob("*.trees"), key=lambda ts: int(ts.stem)) # Iteration order must match inds_file
    for ts in tqdm(ts_lst):
        i = int(ts.stem)
        ts = tskit.load(ts)
        
        if mode == "sub_ts":
            ordered_windows = get_ordered_windows(ts, window_length, A=0, B=1)
            total[key][i] = ordered_windows.tolist()
        
        elif mode == "ts":
            # Recapitation needed, but not mutations!
            ts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-6)
            # ts = msprime.sim_mutations(ts, rate=5e-7)

            total[key][i] = {}
            for j, (A, B) in enumerate(ind_pairs[i]):
                ordered_windows = get_ordered_windows(ts, window_length, A, B)
                total[key][i][j] = ordered_windows.tolist()

with open(out_file, "w") as f:
    json.dump(total, f)
