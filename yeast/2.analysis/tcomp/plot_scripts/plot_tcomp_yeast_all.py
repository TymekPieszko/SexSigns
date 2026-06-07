from multiprocessing import Pool
import os, subprocess, sys, json

# python plot_tcomp_yeast_all.py 1 pair True

PY_SCRIPT = "plot_tcomp_yeast.py"

def run_one(args):
    ploidy, strategy, sc, clade, chrom_id = args
    subprocess.run(
        [
            sys.executable,
            PY_SCRIPT,
            str(ploidy),
            strategy,
            sc,
            clade,
            str(chrom_id),
        ],
        check=True,
    )

def main():
    if len(sys.argv) != 4:
        print("Usage: python plot_tcomp_yeast_all.py 1|2 pair True|False")
        sys.exit(1)

    ploidy = int(sys.argv[1])
    strategy = sys.argv[2]
    sc = sys.argv[3]

    with open("../../../config/clade_dict.json") as f:
        clade_dict = json.load(f)
        clades = list(clade_dict.keys())

    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))

    for chrom_id in range(1, 17):
        args_list = [(ploidy, strategy, sc, clade, chrom_id) for clade in clades]
        with Pool(processes=min(nproc, len(clades))) as pool:
            pool.map(run_one, args_list)

if __name__ == "__main__":
    main()