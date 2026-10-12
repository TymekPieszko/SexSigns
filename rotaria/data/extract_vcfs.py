import shutil
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path
from tqdm import tqdm

####################
### Run as:
# python extract_vcfs.py Rsoc 12

# Fail here if bcftools is not on path
assert shutil.which("bcftools"), "bcftools not on PATH!"

####################
### The worker to parallelise
def extract_vcfs(args):
    bed, vcf, out_dir = args
    pair = bed.stem; s1, s2 = pair.split("~")
    out = out_dir / f"{pair}.vcf.gz"
    # -c1 -C3 ensures only variants sites in the pair get written;
    # invariant sites get redundantly filtered out later anyway.
    subprocess.run(["bcftools", "view", str(vcf), "-R", str(bed), "-s", f"{s1},{s2}",
                    "-c1", "-C3", "-Oz", "-o", str(out)], check=True)
    subprocess.run(["bcftools", "index", "-f", str(out)], check=True)

####################
def main():
    sp = sys.argv[1]
    n_proc = int(sys.argv[2])
    vcf = Path(f"./vcfs_total/{sp}_polar.SNPs.vcf.gz")
    bed_dir = Path(f"./beds_per_pair/{sp}")
    out_dir = Path(f"./vcfs_per_pair/{sp}")

    out_dir.mkdir(parents=True, exist_ok=True)
    beds = sorted(bed_dir.glob("*.bed"))
    print(f"{sp}: {len(beds)} beds -> {out_dir}", flush=True)

    args_list = [(bed, vcf, out_dir) for bed in beds]
    with Pool(n_proc) as pool:
        for _ in tqdm(pool.imap_unordered(extract_vcfs, args_list),
                      total=len(args_list)):
            pass

    print(f"Done: {len(beds)} vcfs -> {out_dir}", flush=True)

if __name__ == "__main__":
    main()
