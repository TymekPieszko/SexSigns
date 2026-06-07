from multiprocessing import Pool
import os, subprocess, sys, json

# python plot_stcomp_yeast_combined_all.py 1 pair True 0.0

PY_SCRIPT = "plot_stcomp_yeast_combined.py"

def run_one(args):
    ploidy, strategy, sc, min_pop_prop, clade = args
    subprocess.run(
        [
            sys.executable,
            PY_SCRIPT,
            str(ploidy),
            strategy,
            sc,
            str(min_pop_prop),
            clade,
        ],
        check=True,
    )

def main():
    if len(sys.argv) != 5:
        print("Usage: python plot_stcomp_yeast_combined_all.py 1|2 pair True|False min_pop_prop")
        sys.exit(1)

    ploidy = int(sys.argv[1])
    strategy = sys.argv[2]
    sc = sys.argv[3]
    min_pop_prop = float(sys.argv[4])

    with open("../../../config/clade_dict.json") as f:
        clade_dict = json.load(f)
        clades = list(clade_dict.keys())

    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))

    args_list = [(ploidy, strategy, sc, min_pop_prop, clade) for clade in clades]

    with Pool(processes=min(nproc, len(clades))) as pool:
        pool.map(run_one, args_list)

if __name__ == "__main__":
    main()