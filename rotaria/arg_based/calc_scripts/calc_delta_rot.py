from pathlib import Path
import numpy as np
from itertools import combinations
from tqdm import tqdm
from multiprocessing import Pool
import tskit, json, sys, allel
from rotspec_utils.calc import polarise_genos, get_info_sites, get_tree_intervals, filter_by_dbl, calc_delta

##############
### Run as:
# python calc_delta_rot.py Rmag 15 10
# python calc_delta_rot.py Rsoc 15 10

##############
### Input:
sp = sys.argv[1]
min_dbl = int(sys.argv[2]) # Min doubletons per CL interval (tree support) as a measure of tree support
n_proc = int(sys.argv[3])
in_dir = Path(f"../ts/{sp}")
vcf_dir = Path(f"../../data/vcfs_per_pair/{sp}")
with open(f"../../data/vcfs_total/{sp}_samples.txt") as f:
    samples = [l.strip() for l in f if l.strip()]
# exclude = ["RS1", "RM9", "RM15"]
exclude = []

##############
### Output:
out_dir = Path(f"../delta/{sp}")
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"d~{min_dbl}.json"
stats_path = out_dir / f"d~{min_dbl}.stats.tsv"
span = lambda d: sum(e - s for ivs in d.values() for s, e in ivs)

##############
# I apologise to myself, or anyone else reading the code, for this solution;
# this basically repeats all the logic from run_sticcs_rot.py,
# but keeps singletons (include_sgl=True), which delta needs.
# I believe this computation could be more elegant if:
# (A) There was a 1-to-1 correspondence between phased segments and VCFs.
# (B) Phased segment genotypes were stored in some custom format.
# I will try to implement (B) in future analyses.

def get_sites(vcf):
    callset = allel.read_vcf(
        str(vcf),
        fields=["variants/CHROM", "variants/POS", "calldata/GT", "calldata/PS",
                "variants/ALT", "variants/AA"]
        )
    contigs = callset["variants/CHROM"]
    pos = np.asarray(callset["variants/POS"])
    genos = allel.GenotypeArray(callset["calldata/GT"])
    # Immediately set phasing mask;
    # unphased sites are marked with 0.
    genos.is_phased = callset["calldata/PS"] > 0
    alt = np.asarray(callset["variants/ALT"])[:,0] # ALT has 3 slots
    anc = callset["variants/AA"]; anc[anc == ""] = "N"

    # Polarisation & info site filter (singletons included: delta needs them)
    genos, mask_polar = polarise_genos(genos, alt, anc)
    pos, contigs = pos[mask_polar], contigs[mask_polar]
    genos, mask_info = get_info_sites(genos, include_sgl=True)
    pos, contigs = pos[mask_info], contigs[mask_info]

    return genos, pos, contigs

##############
### Main:
def worker(pair):
    s1, s2 = pair
    ts_dir = in_dir / f"{s1}~{s2}"
    vcf = vcf_dir / f"{s1}~{s2}.vcf.gz"

    # Once per pair; polarised, phased, biallelic sites incl. singletons
    genos, pos, contigs = get_sites(vcf)
    # Reshape into tskit-style genos
    genos = genos.reshape(-1, 4)

    weighted_sum = 0.0
    total_span = 0.0
    CL_total = 0; CL_retained = 0

    for ts_path in sorted(ts_dir.glob("*.trees")):
        contig, seg, _ = ts_path.stem.rsplit(".", 2)
        s, e = (int(x) for x in seg.split("-"))
        mask_i = (contigs == contig) & (pos >= s+1) & (pos <= e)
        genos_i, pos_i = genos[mask_i], pos[mask_i]-s # Note pos shift
        if genos_i.shape[0] == 0:
            continue

        ts = tskit.load(ts_path)
        interval_dict = get_tree_intervals(ts, include_SX=False)
        CL_total += span(interval_dict)
        interval_dict = filter_by_dbl(interval_dict, genos_i, pos_i, min_dbl)
        CL_retained += span(interval_dict)
        delta_i, span_i = calc_delta(genos_i, pos_i, interval_dict)
        if span_i > 0 and np.isfinite(delta_i):
            weighted_sum += delta_i * span_i
            total_span += span_i

    delta = weighted_sum / total_span if total_span > 0 else np.nan
    print(f"{s1}~{s2} - delta: {delta:.6f}  span: {total_span:.0f}", flush=True)
    return f"{s1}~{s2}", delta, f"{s1}\t{s2}\t{CL_total}\t{CL_retained}\t{100 * CL_retained / CL_total if CL_total else float('nan'):.2f}"

if __name__ == "__main__":
    total = {}
    stats = []
    pairs = [p for p in combinations(samples, 2) if not np.isin(p, exclude).any()]
    with Pool(processes=n_proc) as pool:
        for key, delta, line in tqdm(pool.imap(worker, pairs), total=len(pairs)):
            total[key] = delta
            stats.append(line)

    # Don't write any output if total is empty!
    if len(total) != 0:
        with open(out_path, "w") as f:
            json.dump(total, f, sort_keys=False)
        with open(stats_path, "w") as f:
            f.write("ind1\tind2\tCL_total\tCL_retained\t%_filtered\n" + "\n".join(stats) + "\n")

