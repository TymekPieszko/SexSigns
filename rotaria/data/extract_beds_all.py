from multiprocessing import Pool
from itertools import combinations
from tqdm import tqdm
import subprocess, sys

####################
### Run as:
# python extract_beds_all.py Rmag 3 10
# python -extract_beds_all.py Rsoc 3 9

py_script = "extract_beds.py"

def worker(args):
    sp, s1, s2, min_info = args
    subprocess.run(
        [
            sys.executable,
            "-W", "ignore::UserWarning",
            py_script,
            str(sp),
            str(s1),
            str(s2),
            str(min_info),
        ],
        check=True,
    )

def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_beds_all.py Rmag|Rsoc min_info n_proc")
        sys.exit(1)

    sp = str(sys.argv[1])
    min_info = int(sys.argv[2])
    n_proc = int(sys.argv[3])

    with open(f"./vcfs_total/{sp}_samples.txt") as f:
        samples = [l.strip() for l in f if l.strip()]
    pairs = list(combinations(samples, 2))

    args_list = [(sp, s1, s2, min_info) for s1, s2 in pairs]
    with Pool(processes=n_proc) as pool:
        for _ in tqdm(pool.imap_unordered(worker, args_list),
                      total=len(args_list)):
            pass

if __name__ == "__main__":
    main()