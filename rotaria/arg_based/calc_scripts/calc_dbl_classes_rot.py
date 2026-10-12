from pathlib import Path
import numpy as np
from itertools import combinations
from tqdm import tqdm
import sys, allel
from rotspec_utils.calc import polarise_genos, rm_silwood, rs_silwood

##############
### Run as:
# python calc_dbl_classes_rot.py Rmag silwood
# python calc_dbl_classes_rot.py Rsoc silwood

##############
### Input:
sp = sys.argv[1]
mode = str(sys.argv[2]) # silwood or total
vcf = Path(f"../../data/vcfs_total/{sp}_polar.SNPs.vcf.gz")

##############
### Output:
out_dir = Path(f"../dbl_classes")
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"{sp}_dbl_classes.{mode}.tsv"

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

pairs = list(combinations(samples, 2))
if mode == "silwood":
    silwood = {"Rmag": rm_silwood, "Rsoc": rs_silwood}[sp]
    pairs = [p for p in pairs if set(p) <= set(silwood)]

rows = []
for s1, s2 in tqdm(pairs):
    genos_ij = genos[:, [samples.index(s1), samples.index(s2)]]
    # Exclude sites with missingness or a third allele in the pair
    keep = genos_ij.is_called().all(axis=1) & (genos_ij.max(axis=(1, 2)) <= 1)
    genos_ij, _ = polarise_genos(genos_ij[keep], alt[keep], anc[keep])
    n_der = genos_ij.to_n_alt()
    is_dbl = n_der.sum(axis=1) == 2
    hom = int((is_dbl & (n_der.max(axis=1) == 2)).sum()) # 1/1 + 0/0
    het = int((is_dbl & (n_der.min(axis=1) == 1)).sum()) # 0/1 + 0/1
    rows.append([s1, s2, hom, het, 100 * hom / (hom + het), 100 * het / (hom + het)])

with open(out_path, "w") as f:
    f.write("sample_1\tsample_2\thom\thet\t%_hom\t%_het\n")
    for s1, s2, hom, het, p_hom, p_het in rows:
        f.write(f"{s1}\t{s2}\t{hom}\t{het}\t{p_hom:.3f}\t{p_het:.3f}\n")
    means = np.mean([r[2:] for r in rows], axis=0)
    f.write(f"Mean:\t\t{means[0]:.0f}\t{means[1]:.0f}\t{means[2]:.3f}\t{means[3]:.3f}\n")
