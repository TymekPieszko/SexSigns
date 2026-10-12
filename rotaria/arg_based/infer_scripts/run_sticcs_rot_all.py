from multiprocessing import Pool
from pathlib import Path
from tqdm import tqdm
import subprocess, sys

####################
### Run as:
# python run_sticcs_rot_all.py Rsoc 3 12

py_script = "run_sticcs_rot.py"

def worker(args):
    sp, pair, min_info = args
    subprocess.run(
        [
            sys.executable,
            "-W", "ignore::UserWarning",
            py_script,
            str(sp),
            str(pair),
            str(min_info),
        ],
        check=True,
    )

def main():
    if len(sys.argv) != 4:
        print("Usage: python run_sticcs_rot_all.py Rmag|Rsoc min_info n_proc")
        sys.exit(1)

    sp = str(sys.argv[1])
    min_info = int(sys.argv[2])
    n_proc = int(sys.argv[3])

    vcf_dir = Path(f"../../data/vcfs_per_pair/{sp}")
    pairs = sorted(p.name.removesuffix(".vcf.gz") for p in vcf_dir.glob("*.vcf.gz"))
    print(f"{sp}: {len(pairs)} pairs, min_info={min_info}", flush=True)

    args_list = [(sp, pair, min_info) for pair in pairs]
    with Pool(processes=n_proc) as pool:
        for _ in tqdm(pool.imap_unordered(worker, args_list),
                      total=len(args_list)):
            pass

if __name__ == "__main__":
    main()