from pathlib import Path
import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import groupby
import json, sys
from sexsigns_utils.yeast import get_features

# python get_callable_regions.py 12.West_African_cocoa
# python get_callable_regions.py '1.Wine_European_(subclade_3)'
# python get_callable_regions.py 3.Brazilian_bioethanol
# python get_callable_regions.py 25.Sake

##############
### Input:
clade = sys.argv[1]
clade_dict = f"./clade_dict.json"
clade_dict = json.load(open(clade_dict))
features = "../config/chrom_features.tsv" # 
features = pd.read_csv(features, sep="\t")
align_dir = Path(f"/data/biol-bdelloids/scro4331/SexSigns/yeast/0.calls/PH/0.1.alignments/")

##############
### Output:
out = Path(f"./callable_regions/{clade}.tsv")
out.parent.mkdir(exist_ok=True)
out = open(out, "w")
out.write("sample\tchrom_id\tstart\tend\n")

def get_coverage(paf, chrom_len):
    # chrom_len is used to safeguard against empty pafs
    coverage = None
    with open(paf) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.split("\t")

            if coverage is None:
                chrom_len_from_paf = int(fields[6])
                assert chrom_len_from_paf == chrom_len # Just to be sure
                coverage = [0] * chrom_len_from_paf

            start = int(fields[7])
            end   = int(fields[8])

            for pos in range(start, end):
                coverage[pos] += 1
    if coverage is None:
        return np.array([0] * chrom_len)        
    return np.array(coverage)

def get_callable_regions(cov_per_pos, target_cov):
    out = []
    i = 0
    for k,v in groupby(cov_per_pos):
        increment = len(list(v))
        interval = (i, i+increment)
        i += increment
        if k != target_cov:
            continue
        out.append(interval)
    return np.array(out, dtype=int)

def subtract_interval(intervals, to_subtract):
    s2, e2 = to_subtract
    out = []
    for s1, e1 in intervals:
        if e1 <= s2 or s1 >= e2:
            out.append((s1, e1))
        else:
            if s1 < s2:
                out.append((s1, s2))
            if e1 > e2:
                out.append((e2, e1))
    return np.array(out, dtype=int).reshape(-1, 2) # Safeguard against empty case


for sample in tqdm(clade_dict[clade]):
    for chrom_id in range(1,17):
        chrom_dir = align_dir / f"chrom_{chrom_id}"
        # Get chromosome features
        centr_start, centr_end, chrom_len = get_features(features, chrom_id)
        # Get coverage from alignments
        paf_1 = chrom_dir / f"{sample}.HP1.paf"
        paf_2 = chrom_dir / f"{sample}.HP2.paf"
        cov_1 = get_coverage(paf_1, chrom_len)
        cov_2 = get_coverage(paf_2, chrom_len)
        # Combine coverage from haplotypes 1 and 2
        cov = cov_1 + cov_2
        callable_regions = get_callable_regions(cov, target_cov=2)
        # Subtract the centromere region
        callable_regions = subtract_interval(callable_regions, [centr_start, centr_end])
        for s, e in callable_regions:
            # Skip regions < 1 kb
            if e-s < 1000:
                continue
            out.write(f"{sample}\t{chrom_id}\t{s}\t{e}\n")
out.close()