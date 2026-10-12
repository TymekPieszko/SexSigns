from pathlib import Path
import numpy as np
from itertools import combinations
from multiprocessing import Pool
import tskit, sys
from sexsigns_utils.calc import subsample_to_inds, get_biallelic_01, get_tree_intervals, calc_delta, filter_by_dbl, calc_kappa, codes, calc_tcomp_vRanks, calc_category_props

##############
### Run as:
# python calc_stats_for_abc.py GC 10 0.1 15 16
# python calc_stats_for_abc.py GC 8 0.1 15 16

# Note that INDS~8 uses the first 8 individuals of the INDS~10 outputs.

##############
### Input:
model = sys.argv[1]
n_inds = int(sys.argv[2])
bias = float(sys.argv[3])
min_dbl = int(sys.argv[4])
n_proc = int(sys.argv[5])
window_size = 5000
src_inds = 2 if n_inds == 2 else 10
sub_dir = Path(f"../sim_pipeline/sim_output/{model}/1.0.sub_ts/INDS~{src_inds}")
sticcs_dir = Path(f"../sim_pipeline/sim_output/{model}/2.0.sticcs/INDS~{n_inds}")
pairs = list(combinations(range(n_inds), 2))

##############
### Output:
out_file = Path(f"./for_abc/{model}_stats_inds~{n_inds}_b~{bias}_d~{min_dbl}.tsv")
out_file.parent.mkdir(exist_ok=True)

##############
def get_biased_windows(ts, window_size, bias):
    windows = np.arange(0, ts.sequence_length + window_size, window_size)
    het = ts.diversity([[0, 1], [2, 3]], windows=windows, mode="site").mean(axis=1)
    top = np.sort(np.argsort(het)[-int(len(het) * bias):])
    return np.column_stack([windows[:-1], windows[1:]])[top]

def calc_pair_stats(ts, ts_sticcs):
    genos, pos = get_biallelic_01(ts.genotype_matrix(), ts.sites_position)
    kappa = calc_kappa(genos)
    windows = get_biased_windows(ts, window_size, bias)
    ts_sticcs = ts_sticcs.keep_intervals(windows, simplify=False)
    props = calc_category_props(calc_tcomp_vRanks(ts_sticcs, {code: 0.0 for code in codes}))
    delta, _ = calc_delta(genos, pos, filter_by_dbl(get_tree_intervals(ts_sticcs, include_SX=False), genos, pos, min_dbl))
    return kappa, props["CL"], props["AR"], props["SX"], delta

def worker(scen_dir):
    sex, rec = (p.split("~")[1] for p in scen_dir.parts[-3:-1])
    rows = []
    for ts_path in sorted(scen_dir.glob("*.trees"), key=lambda p: int(p.stem)):
        rep_dir = sticcs_dir / f"SEX~{sex}" / f"REC~{rec}" / scen_dir.name / ts_path.stem
        if not all((rep_dir / f"{i}~{j}.trees").exists() for i, j in pairs):
            continue
        ts = tskit.load(ts_path)
        stats = [calc_pair_stats(subsample_to_inds(ts, i, j), tskit.load(rep_dir / f"{i}~{j}.trees")) for i, j in pairs]
        rows.append(f"{sex}\t{rec}\t" + "\t".join(f"{v:.6f}" for v in np.nanmean(stats, axis=0)))
    print(sex, rec, len(rows), flush=True)
    return rows

##############
### Main:
if __name__ == "__main__":
    scen_dirs = sorted(sub_dir.glob("SEX~*/REC~*/MUT~*"))
    with Pool(processes=n_proc) as pool:
        rows = [row for rows in pool.imap_unordered(worker, scen_dirs) for row in rows]
    with open(out_file, "w") as f:
        f.write("sex\trec\tkappa\tCL\tAR\tSX\tdelta\n" + "\n".join(rows) + "\n")
