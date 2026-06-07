from pathlib import Path
import numpy as np
from tqdm import tqdm
import allel, json, sys, os
from multiprocessing import Pool
from sexsigns_utils.calc import to_tskit, to_allel, get_biallelic_01
from sexsigns_utils.yeast import get_callable_regions, get_callable_mask

# python calc_ld.py 12.West_African_cocoa 1000 20000 100 100

### !!! THIS SCRIPT NEEDS TO BE ORGANISED !!! ###

##############
### Input:
clade = sys.argv[1]
k = int(sys.argv[2])
d_max = int(sys.argv[3])
window_size = int(sys.argv[4])
reps = int(sys.argv[5])

vcf_base = Path("/data/biol-bdelloids/scro4331/SexSigns/yeast/0.calls/PH/0.4.total_vcfs")
clade_dict = "../../../config/clade_dict.json"

callable_df = Path(f"../../../config/callable_regions/{clade}.tsv")

##############
### Output:
out_dir = Path("../LD") / clade
out_dir.mkdir(exist_ok=True)

##############
# Get target samples
with open(clade_dict) as f:
    clade_dict = json.load(f)

targets = clade_dict[clade]
print(targets)

##############
# Load callable regions
import pandas as pd
callable_df = pd.read_csv(callable_df, sep="\t")

##############
# Load ALL chromosomes once

all_genos = []
all_pos = []
offset = 0

for chrom in range(1, 17):

    vcf = vcf_base / f"chrom_{chrom}" / f"chrom_{chrom}.biSNPs.polar.vcf.gz"

    callset = allel.read_vcf(
        str(vcf),
        samples=targets,
        fields=['variants/POS', 'calldata/GT']
    )

    genos = allel.GenotypeArray(callset['calldata/GT'])
    pos = callset['variants/POS']

    genos, pos = get_biallelic_01(to_tskit(genos), pos)
    genos = to_allel(genos)

    #######################
    ### Apply callable mask
    callable_regions_per_sample = [
        get_callable_regions(callable_df, chrom, target)
        for target in targets
    ]

    callable_mask_per_sample = get_callable_mask(
        callable_regions_per_sample,
        pos
    )

    genos.mask = callable_mask_per_sample
    #######################

    pos = pos + offset

    all_genos.append(genos)
    all_pos.append(pos)

    offset += pos.max() + 1

genos = allel.GenotypeArray(np.concatenate(all_genos, axis=0))
pos = np.concatenate(all_pos)

##############

def calc_r2(genos, pos, k, d_max):
    r2_total = []
    dist_total = []

    for _ in range(k):

        i = np.random.randint(0, len(genos) - 1)
        g_i = genos[i]
        p_i = pos[i]

        g_other = allel.GenotypeArray(np.delete(genos, i, axis=0))
        p_other = np.delete(pos, i)

        alt_count_i = [g_i.to_n_alt()]
        alt_counts = g_other.to_n_alt()

        r = np.atleast_1d(
            allel.rogers_huff_r_between(alt_count_i, alt_counts)
        )[0]

        r2 = r ** 2
        dist = np.abs(p_other - p_i)

        mask = (dist != 0) & (dist < d_max)

        r2_total.extend(r2[mask])
        dist_total.extend(dist[mask])

    r2_total = np.array(r2_total)
    dist_total = np.array(dist_total)

    sort_idx = np.argsort(dist_total)
    return r2_total[sort_idx], dist_total[sort_idx]

##############
# Parallel replicate function

def run_rep(rep):

    r2, dist = calc_r2(genos, pos, k, d_max)

    r2_windows = allel.windowed_statistic(
        dist,
        r2,
        statistic=lambda x: np.nanmean(x),
        start=0,
        stop=int(dist[-1]),
        size=window_size
    )[0]

    out_file = out_dir / f"{rep}.txt"
    np.savetxt(out_file, r2_windows, fmt="%.3f")


##############

if __name__ == "__main__":

    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
    print(f"Using {nproc} processes (parallel over replicates)")

    with Pool(processes=min(nproc, reps)) as pool:
        list(tqdm(pool.imap_unordered(run_rep, range(reps)), total=reps))