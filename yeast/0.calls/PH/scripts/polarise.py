from cyvcf2 import VCF, Writer
from pathlib import Path
import numpy as np
import sys, tqdm

"""
Logic:
    - SNP in vcf_in BUT NOT vcf_outgroup - site where outgroup samples are fixed for the reference allele. Nothing to be done.
    - SNP shared between vcf_in and vcf_outgroup:
        - ALT count >= 5: ALT assumed ancestral, flip REF/ALT and invert genotypes
        - ALT count 2-4: ambiguous, exclude variant
        - ALT count <= 1: REF assumed ancestral, keep as is
"""

###############
### Inputs:
# Outgroup VCFs were obtained as described in /data/biol-bdelloids/scro4331/SexSigns_2025/yeast/0.calls/PH/ancestral_states/commands.sh
vcf_in = sys.argv[1]
vcf_outgroup = sys.argv[2]
# Set of samples to polarise against;
# chose the 3 samples of the Taiwanese clade (clade 17),
# as defined by Peter et al. (2018) (https://doi.org/10.1038/s41586-018-0030-5).
outgroup_samples = ['AMH', 'CEG', 'CEI']

###############
### Outputs:
vcf_out = Path(sys.argv[3])
log_file = vcf_out.parent / f"{vcf_out.name}.log" 

###############
### Main:

### Build a look-up dict for test VCF
vcf_outgroup = VCF(vcf_outgroup)
outgroup_idx = [vcf_outgroup.samples.index(s) for s in outgroup_samples]
test_dict = {}
for v in vcf_outgroup:
    pos = v.POS
    ref = v.REF
    alt = v.ALT[0]

    outgroup_genos = np.array([v.genotypes[i][:2] for i in outgroup_idx])
    alt_count = np.sum(outgroup_genos == 1)

    test_dict[pos] = [ref, alt, alt_count]

### Polarise the target
vcf_in = VCF(vcf_in)
out = Writer(vcf_out, vcf_in)
total = 0
ref_out_fixed = 0
matched = 0
flipped = 0
correct = 0
excluded = 0
for v in tqdm.tqdm(vcf_in):
    # Sanity filter for biallelic SNPs.
    if not v.is_snp or len(v.ALT) != 1:
        continue
    total += 1
    
    p = v.POS
    r = v.REF
    a = v.ALT[0]
    if (p in test_dict) and (set([r, a]) == set(test_dict[p][:2])):
        matched += 1
        # Consider the number of alternative alleles in outgroup
        alt_count = test_dict[p][2]
        if alt_count >= 5:
            # Flip polarity!
            flipped += 1
            #######################
            v.REF, v.ALT = v.ALT[0], [v.REF]
            for g in v.genotypes:
                for i in (0, 1):
                    if g[i] == 0:
                        g[i] = 1
                    elif g[i] == 1:
                        g[i] = 0
            #######################
        elif alt_count <= 1:
            # Polarity is correct!
            correct += 1
        else:
            # Exclude variant!
            excluded += 1
            continue
    else:
        ref_out_fixed += 1
        pass
    v.genotypes = v.genotypes
    out.write_record(v)
out.close()

# print(f"Total SNPs in vcf_outgroup: {len(test_dict)}")
# print(f"Total SNPs in vcf_in: {total}")
# print(f"Ref + outgroup fixed for the same allele: {ref_out_fixed}")
# print("Now, for SNPs in target with matches in test...")
# print(f"SNPs with matching pos + alleles in vcf_PH: {matched}")
# print(f"Out of which {flipped} were flipped, {correct} were correct and {excluded} were excluded.")

with open(log_file, "w") as log:
    log.write(f"Total SNPs in vcf_outgroup: {len(test_dict)}\n")
    log.write(f"Total SNPs in vcf_in: {total}\n")
    log.write(f"Ref + outgroup fixed for the same allele: {ref_out_fixed}\n")
    log.write("Now, for SNPs in target with matches in test...\n")
    log.write(f"SNPs with matching pos + alleles in vcf_PH: {matched}\n")
    log.write(f"Out of which {flipped} were flipped, {correct} were correct and {excluded} were excluded.\n")