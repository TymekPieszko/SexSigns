import numpy as np
import allel
import sys
import os
from multiprocessing import Pool
from tqdm import tqdm
import unphased_local as uph

# python calc_ld.py ../../data/vcfs_total/Rmag_polar.SNPs.vcf.gz rm_SIL 10000 10000 100 100 10
# python calc_ld.py ../../data/vcfs_total/Rmag_polar.SNPs.vcf.gz rm_total 10000 10000 100 100 10
# python calc_ld.py ../../data/vcfs_total/Rsoc_polar.SNPs.vcf.gz rs_SIL 10000 10000 100 100 10
# python calc_ld.py ../../data/vcfs_total/Rsoc_polar.SNPs.vcf.gz rs_total 10000 10000 100 100 10
vcf = sys.argv[1]
pop = str(sys.argv[2])
K = int(sys.argv[3])
d_max = int(sys.argv[4])
window_size = int(sys.argv[5]) # i.e., d_max / window_size windows
n_reps = int(sys.argv[6]) # bootstrap replicates (one output file each)
n_proc = int(sys.argv[7])

if pop == "rm_SIL":
    samples = uph.rm_SIL_samples
if pop == "rm_total":
    samples = uph.rm_SIL_samples + uph.rm_other_samples
if pop == "rs_SIL":
    samples = uph.rs_SIL_samples
if pop == "rs_total":
    samples = uph.rs_SIL_samples + uph.rs_other_samples
out_dir = "../ld/" + pop + "_K_" + f"{str(K)}" + "_window_" + str(window_size)
callset = allel.read_vcf(vcf, samples=samples, fields=['variants/CHROM', 'variants/POS', 'calldata/GT'])
genos = allel.GenotypeArray(callset['calldata/GT'])
pos =  callset['variants/POS']
scaffolds = callset['variants/CHROM']

def calc_r2(genos, pos, scaffolds, K, d_max):
    # K = number of target SNPs
    # d_max = maximum distance from target SNP
    r2_total = []
    dist_total = []
    for _ in range(K):
        ##############
        target = np.random.randint(0, len(genos)-1)
        target_geno = genos[target]
        target_pos = pos[target]
        target_scaffold = scaffolds[target]
        ##############
        other = np.where(scaffolds == target_scaffold)[0]
        other_genos = genos[other]
        other_pos = pos[other]
        # Skip if target is the only SNP at a scaffold
        if len(other) == 1:
            continue
        ##############
        target_geno_n = [target_geno.to_n_alt()]
        other_genos_n = other_genos.to_n_alt()
        r_i = np.atleast_1d(allel.rogers_huff_r_between(target_geno_n, other_genos_n))[0]
        r2_i = r_i ** 2
        dist_i = np.abs(other_pos - target_pos)
        mask = (dist_i!=0) & (dist_i<d_max)
        r2_total.extend(r2_i[mask])
        dist_total.extend(dist_i[mask])
    # Sort by distance
    r2_total = np.array(r2_total)
    dist_total = np.array(dist_total)
    sort_idx = np.argsort(dist_total)
    return r2_total[sort_idx], dist_total[sort_idx]

def worker(i):
    np.random.seed(i) # replicates must differ even if workers are forked with a shared RNG state
    r2, dist = calc_r2(genos, pos, scaffolds, K, d_max)
    r2_windows = allel.windowed_statistic(dist, r2, statistic=lambda x : np.nanmean(x), start=0, stop=d_max, size=window_size)[0]
    np.savetxt(out_dir + f"/{i}.txt", r2_windows, fmt="%.3f")

if __name__ == "__main__":
    os.makedirs(out_dir, exist_ok=True)
    with Pool(processes=n_proc) as pool:
        list(tqdm(pool.imap_unordered(worker, range(n_reps)), total=n_reps))
