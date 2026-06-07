from pathlib import Path
from tqdm import tqdm
import tskit, json, sys
from sexsigns_utils.calc import codes, calc_tcomp_vTracked
from sexsigns_utils.yeast import calc_stcomp
import numpy as np

##############
### About:
# Calculate tree composition
# Run for all clades / chromosomes using calc_tcomp_yeast_all.py
# Tymek Pieszko; 2026-04-14

##############
### Run as:
# python calc_tcomp_yeast.py 1 pair True 12.West_African_cocoa 1
# python calc_tcomp_yeast.py 1 pair True '1.Wine_European_(subclade_3)' 4

##############
### Input:
ploidy = int(sys.argv[1])
strategy = sys.argv[2]
sc = sys.argv[3]
clade = sys.argv[4]
chrom_id = int(sys.argv[5])
ts_dir = Path(f"../../../1.args/sticcs/ploidy~{ploidy}_strategy~{strategy}_sc~{sc}/chrom_{chrom_id}/{clade}")
##############
### Output:
out_dir = Path(f"../data/ploidy~{ploidy}_strategy~{strategy}_sc~{sc}/chrom_{chrom_id}/")
out_dir.mkdir(exist_ok=True, parents=True)
out_file_1 = f"{clade}.json"
out_file_2 = f"{clade}.spatial.json"
out_path_1 = out_dir / out_file_1
out_path_2 = out_dir / out_file_2

total_1 = {}
total_2 = {}
for ts in tqdm(ts_dir.glob("*.trees")):
    pair = str(ts.stem)
    ts = tskit.load(ts)

    # Prepare dicts to store tcomp & (spatial) stcomp
    tcomp = {code : 0.0 for code in codes}
    stcomp = {code : [] for code in codes + ["poly", "multi_root"]}

    # Calculate tcomp & stcomp
    tcomp = calc_tcomp_vTracked(ts, tcomp)
    stcomp = calc_stcomp(ts, stcomp)

    total_1[pair] = tcomp
    total_2[pair] = stcomp

# Don't write any output if total is empty!
if len(total_1) != 0:
    with open(out_path_1, 'w') as f:
        json.dump(total_1, f, sort_keys=False)
if len(total_2) != 0:
    with open(out_path_2, 'w') as f:
        json.dump(total_2, f, sort_keys=False)