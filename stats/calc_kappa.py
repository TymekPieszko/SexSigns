from pathlib import Path
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
import tskit, msprime, pyslim, json, sys, warnings
from sexsigns_utils.calc import read_inds, get_biallelic_01, calc_kappa

warnings.simplefilter('ignore')

#############
### ABOUT ###
#############
# Calculate kappa for target individuals drawn from simulation replicates. 

### Run as:
# python calc_kappa.py GC 16

###################
### Input:
sim_model = sys.argv[1]
n_workers = int(sys.argv[2])
inds_file = "./inds_60.tsv"
target = Path(f"../sim_pipeline_2/sim_output/{sim_model}/0.ts/")

###################
### Output:
out_file = Path(f"./kappa/{sim_model}.json")
out_file.parent.mkdir(exist_ok=True)

###################
def worker(args):
    ts_file, ind_pairs = args
    i = int(ts_file.stem)

    ts = tskit.load(ts_file)
    ts = pyslim.recapitate(ts, ancestral_Ne=1000, recombination_rate=1e-6)
    ts = msprime.sim_mutations(ts, rate=5e-7)
    genos = ts.genotype_matrix()
    pos = ts.sites_position
    genos = get_biallelic_01(genos, pos)[0]

    kappa_j_lst = []
    for ind_1, ind_2 in ind_pairs:
        nodes_1 = ts.individual(ind_1).nodes
        nodes_2 = ts.individual(ind_2).nodes
        targets = list(nodes_1) + list(nodes_2)

        genos_AB = genos[:, targets]
        kappa_j = calc_kappa(genos_AB)
        kappa_j_lst.append(kappa_j)

    return i, np.nanmean(kappa_j_lst)

###################
### Main:

ind_pairs = read_inds(inds_file)
total = {}
for dir in target.glob("SEX~*/REC~*"):
    sex = float(dir.parts[-2].split("~")[1])
    rec = float(dir.parts[-1].split("~")[1])
    key = f"{sex}_{rec}"
    print(key)

    ts_lst = sorted(dir.glob("*.trees"), key=lambda f: int(f.stem))
    tasks = [(f, ind_pairs[int(f.stem)]) for f in ts_lst]
    
    kappa_lst = []
    with Pool(processes=n_workers) as pool:
        for i, kappa in tqdm(
            pool.imap_unordered(worker, tasks),
            total=len(tasks),
        ):
            kappa_lst.append(kappa)

    print(np.nanmean(kappa_lst))
    total[key] = kappa_lst

with open(out_file, "w") as f:
    json.dump(total, f)
