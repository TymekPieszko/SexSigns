from pathlib import Path
import numpy as np
from tqdm import tqdm
from rotspec_utils.sticcs import run_sticcs
from rotspec_utils.calc import read_bed, get_info_sites, polarise_genos
import sys, allel

####################
### Run as:
# python run_sticcs_rot.py Rsoc AK27~RSOC406 3

####################
### Inputs:
sp = sys.argv[1]
pair = sys.argv[2]
# Info site filter needs to be repeated from extract_beds.py,
# since the filtered arrays are what sticcs is given.
# BUT IT IS NOT INTENDED to be adjusted here;
# it should be adjusted when calculating tcomp.
min_info = int(sys.argv[3])
vcf = f"../../data/vcfs_per_pair/{sp}/{pair}.vcf.gz"
bed = f"../../data/beds_per_pair/{sp}/{pair}.bed"

####################
### Outputs:
out_dir = Path(f"../ts/{sp}/{pair}")
out_dir.mkdir(parents=True, exist_ok=True)

####################
### Parse VCF
callset = allel.read_vcf(
    str(vcf),
    fields=["samples", "variants/CHROM", "variants/POS", "calldata/GT", "calldata/PS",
            "variants/REF", "variants/ALT", "variants/AA"]
    )
samples = callset["samples"]
contigs = callset["variants/CHROM"]
pos = np.asarray(callset["variants/POS"])
genos = allel.GenotypeArray(callset["calldata/GT"])
# Immediately set phasing mask;
# unphased sites are marked with 0.
genos.is_phased = callset["calldata/PS"] > 0
ref = np.asarray(callset["variants/REF"])
alt = np.asarray(callset["variants/ALT"])[:,0] # ALT has 3 slots
anc = callset["variants/AA"]; anc[anc == ""] = "N"
# print(anc)

####################
### Parse BED
block_info = read_bed(bed)

####################
### Polarasition & info site filter
genos, mask_polar = polarise_genos(genos, alt, anc)
pos, contigs = pos[mask_polar], contigs[mask_polar]
genos, mask_info = get_info_sites(genos)
pos, contigs = pos[mask_info], contigs[mask_info]

####################
### Run sticcs
# Reshape into tskit-style genos
genos = genos.reshape(-1, 4)
for block in block_info:
    contig = block["contig"]
    s, e = int(block["start"]), int(block["end"])
    # Subset to block
    mask_i = (contigs == contig) & (pos >= s+1) & (pos <= e); n_info = mask_i.sum()
    # Enforce minimum number of info sites
    if n_info < min_info:
        continue
    genos_i, pos_i = genos[mask_i], pos[mask_i]-s # Note pos shift
    ts = run_sticcs(
        genos_i,
        pos_i,
        e-s,
        forced_ploidy=1,
        second_chances=True,
        silent=True
    )
    ts.dump(out_dir / f"{contig}.{s}-{e}.n_info~{n_info}.trees")