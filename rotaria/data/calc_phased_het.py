from pathlib import Path
import numpy as np
import sys, allel
from tqdm import tqdm
from rotspec_utils.calc import get_biallelic_01_allel, rm_silwood, rs_silwood

### Run as:
# python calc_phased_het.py silwood 10
# python calc_phased_het.py total 10

###################
### Input:
mode = str(sys.argv[1]) # silwood or total
min_info = int(sys.argv[2])
sp_lst = ["Rmag", "Rsoc"]

###################
### Output:
out_file = Path(f"./phased_het_i~{min_info}.{mode}.tsv")
out_file.parent.mkdir(exist_ok=True)

###################
# This is
def get_segments(pair_dir, min_info):
    segments = {}
    for f in pair_dir.glob("*.trees"):
        contig, coords, n_info = f.stem.rsplit(".", 2)
        if int(n_info.split("~")[1]) < min_info:
            continue
        s, e = map(int, coords.split("-"))
        segments.setdefault(contig, []).append([s + 1, e])
    return {contig: np.array(sorted(ivs)) for contig, ivs in segments.items()}

def calc_phased_het(genos, pos, scaffolds, segments):
    hets = np.zeros(genos.n_samples)
    for contig, windows in segments.items():
        mask = scaffolds == contig
        genos_c, pos_c = get_biallelic_01_allel(genos[mask], pos[mask])
        for i in range(genos.n_samples):
            hets[i] += allel.windowed_statistic(pos_c, genos_c[:, [i]], statistic=lambda x: x.is_het().sum(), windows=windows)[0].sum()
    span = sum((windows[:, 1] - windows[:, 0] + 1).sum() for windows in segments.values())
    return hets / span

###################
### Main:
sp_means = {}
for sp in sp_lst:
    het_lst = []
    silwood = {"Rmag": rm_silwood, "Rsoc": rs_silwood}[sp]
    for pair_dir in tqdm(sorted(Path(f"../arg_based/ts/{sp}").glob("*~*"))):
        pair = pair_dir.name
        if mode == "silwood" and not set(pair.split("~")) <= set(silwood):
            continue
        segments = get_segments(pair_dir, min_info)
        callset = allel.read_vcf(f"./vcfs_per_pair/{sp}/{pair}.vcf.gz", fields=["variants/CHROM", "variants/POS", "calldata/GT"])
        genos = allel.GenotypeArray(callset["calldata/GT"])
        het = calc_phased_het(genos, callset["variants/POS"], callset["variants/CHROM"], segments)
        het_lst.append(het.mean())
    sp_means[sp] = np.mean(het_lst)

with open(out_file, "w") as f:
    for sp, het in sp_means.items():
        f.write(f"{sp}\t{het}\n")
    f.write(f"mean\t{np.mean(list(sp_means.values()))}\n")
