from pathlib import Path
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
import tskit, msprime, pyslim, sys, json, warnings
from sexsigns_utils.calc import read_inds, get_biallelic_01, to_allel, calc_hi, calc_fis, calc_ld_tskit

warnings.filterwarnings("ignore")

#############
### ABOUT ###
#############
# Calculate classical statistics for target individuals drawn from simulation replicates. 

### Run as:
# python calc_classical_stats.py GC hi 16
# python calc_classical_stats.py GC fis 16
# python calc_classical_stats.py GC ld 16
# python calc_classical_stats.py CO hi 16
# python calc_classical_stats.py CO fis 16
# python calc_classical_stats.py CO ld 16

###################
### Input:
model = sys.argv[1]
stat = sys.argv[2]
n_workers = int(sys.argv[3])
inds_file = "./inds_60.tsv"
target = Path(f"../sim_pipeline_2/sim_output/{model}/0.ts/")

###################
### Output:
out_file = f"./{stat}/{model}.json"
Path(out_file).parent.mkdir(exist_ok=True)

###################
def worker(args):
    ts_file, stat, inds = args

    ts = tskit.load(str(ts_file))
    ts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-6)
    nodes = np.concatenate([ts.individual(ind).nodes for pair in inds for ind in pair]) # inds are a list of tuples
    ts = ts.simplify(samples=nodes) # Stats are calculated for all inds at once
    ts = msprime.sim_mutations(ts, rate=5e-7)
    genos = ts.genotype_matrix()
    pos = ts.sites_position
    seq_len = ts.sequence_length
    genos, pos = get_biallelic_01(genos, pos)

    if stat == "hi":
        val = calc_hi(to_allel(genos), seq_len)
    elif stat == "fis":
        val = calc_fis(to_allel(genos))
    # elif stat == "LD_0-1000":
    elif stat == "ld":
        # val = calc_ld(to_allel(genos), pos, 500, 1000) # Takes genos, pos, k, d_max
        val = calc_ld_tskit(ts, 0, 1000) # Takes ts, d1, d2
        # val = calc_fold_ld_decay(to_allel(genos), pos, seq_len, k=1000, d=10000, w=1000)
    return val
###################

###################
### Main:

ind_pairs = read_inds(inds_file)

total_reps = {}
for ts_dir in target.glob("SEX~*/REC~*"):
    sex = float(ts_dir.parts[-2].split("~")[1])
    rec = float(ts_dir.parts[-1].split("~")[1])
    key = f"{sex}_{rec}"
    print(key)

    ts_files = sorted(ts_dir.glob("*.trees"), key=lambda f: float(f.stem))
    tasks = [(f, stat, ind_pairs[int(f.stem)]) for f in ts_files]

    reps = []
    with Pool(processes=n_workers) as pool:
        for val in tqdm(
            pool.imap_unordered(worker, tasks),
            total=len(tasks),
        ):
            reps.append(val)

    total_reps[key] = reps
    print(np.nanmean(reps))

with open(out_file, "w") as f:
    json.dump(total_reps, f)
