from pathlib import Path
from multiprocessing import Pool
import sys, os, json, subprocess

os.environ["PYTHONWARNINGS"] = "ignore"

PY_SCRIPT = "calc_stats.py"

def run_clade(clade):
    subprocess.run(
        [sys.executable, PY_SCRIPT, clade],
        check=True
    )


def main():

    clades_file = Path("./delta_clades.txt")

    with open("../../../config/clade_dict.json") as f:
        clade_dict = json.load(f)
        clades = list(clade_dict.keys())

    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))

    print(f"\nUsing {nproc} processes (parallel over clades)...\n", flush=True)

    with Pool(processes=min(nproc, len(clades))) as pool:
        pool.map(run_clade, clades)


if __name__ == "__main__":
    main()