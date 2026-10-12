from pathlib import Path
import numpy as np
from itertools import combinations
from tqdm import tqdm
import json, sys, allel
from rotspec_utils.calc import polarise_genos, calc_kappa

##############
### Run as:
# python calc_kappa_rot.py Rmag
# python calc_kappa_rot.py Rsoc

##############
### Input:
sp = sys.argv[1]
vcf = Path(f"../../data/vcfs_total/{sp}_polar.SNPs.vcf.gz")

##############
### Output:
out_dir = Path(f"../kappa")
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"{sp}.json"

##############
### Main:
callset = allel.read_vcf(
    str(vcf),
    fields=["samples", "calldata/GT", "variants/ALT", "variants/AA"]
    )
samples = list(callset["samples"])
genos = allel.GenotypeArray(callset["calldata/GT"])
alt = np.asarray(callset["variants/ALT"])[:,0] # ALT has 3 slots
anc = callset["variants/AA"]; anc[anc == ""] = "N"

total = {}
for s1, s2 in tqdm(list(combinations(samples, 2))):
    genos_ij = genos[:, [samples.index(s1), samples.index(s2)]]
    # Exclude sites with missingness or a third allele in the pair
    keep = genos_ij.is_called().all(axis=1) & (genos_ij.max(axis=(1, 2)) <= 1)
    genos_ij, _ = polarise_genos(genos_ij[keep], alt[keep], anc[keep])
    total[f"{s1}~{s2}"] = calc_kappa(genos_ij.reshape(-1, 4))

# Don't write any output if total is empty!
if len(total) != 0:
    with open(out_path, "w") as f:
        json.dump(total, f, sort_keys=False)
