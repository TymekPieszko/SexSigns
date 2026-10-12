from pathlib import Path
import numpy as np
import sys, allel
from rotspec_utils.calc import rm_silwood, rs_silwood

####################
### Run as:
# python calc_snp_stats.py silwood
# python calc_snp_stats.py total

####################
### Inputs:
mode = str(sys.argv[1]) # silwood or total
sp_lst = ["Rmag", "Rsoc"]

####################
### Outputs:
out_file = Path(f"./snp_stats.{mode}.tsv")

####################
### SNPs biallelic and polymorphic among the selected individuals;
### phasing is counted over their heterozygous genotypes.
rows = []
for sp in sp_lst:
    samples = {"Rmag": rm_silwood, "Rsoc": rs_silwood}[sp] if mode == "silwood" else None
    callset = allel.read_vcf(f"./vcfs_total/{sp}_polar.SNPs.vcf.gz", samples=samples,
                             fields=["calldata/GT", "calldata/PS", "variants/AA"])
    genos = allel.GenotypeArray(callset["calldata/GT"])
    keep = genos.count_alleles().is_biallelic_01()
    genos, ps, anc = genos[keep], callset["calldata/PS"][keep], callset["variants/AA"][keep]
    is_polar = np.isin(anc, list("ACGT"))
    is_het = genos.is_het()
    is_phased = is_het & (ps > 0)
    n_het = is_het.sum()
    rows.append([sp, len(genos), n_het, 100 * is_polar.mean(),
                 100 * is_phased.sum() / n_het, 100 * (is_phased & is_polar[:, None]).sum() / n_het])

with open(out_file, "w") as f:
    f.write("species\tn_SNPs\tn_het_genotypes\t%_SNPs_polarised\t%_het_phased\t%_het_phased_polarised\n")
    for sp, n, n_het, *pcts in rows:
        f.write("\t".join([sp, str(n), str(n_het)] + [f"{p:.3f}" for p in pcts]) + "\n")
