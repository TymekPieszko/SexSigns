from pathlib import Path
from tqdm import tqdm
from itertools import combinations
import tskit, json, sys
from rotspec_utils.calc import calc_tcomp_vRanks, codes, calc_category_props

##############
### Run as:
# python calc_tcomp_rot.py Rmag 3
# python calc_tcomp_rot.py Rsoc 3

# Min number of informative sites per segment is recorded in the output directory name.

##############
### Input:
sp = sys.argv[1]
min_info = int(sys.argv[2])
in_dir = Path(f"../ts/{sp}")
with open(f"../../data/vcfs_total/{sp}_samples.txt") as f:
    samples = [l.strip() for l in f if l.strip()]

##############
### Output:
out_dir = Path(f"../tcomp/{sp}/i~{min_info}")
out_dir.mkdir(parents=True, exist_ok=True)

##############
### Main:
pairs = list(combinations(samples, 2))
for s1, s2 in pairs:
    ts_dir = in_dir / f"{s1}~{s2}"
    # print(ts_dir)
    tcomp = {code:0.0 for code in codes}
    for ts_path in tqdm(ts_dir.glob("*.trees")):
        # print(ts_path)
        # print(ts_path.name)
        n_info = int(ts_path.stem.split(".")[-1].split("~")[1])
        # print(n_info)
        if n_info < min_info:
            continue
        ts = tskit.load(ts_path)
        tcomp = calc_tcomp_vRanks(ts, tcomp)

    print(f"{s1}~{s2} - ", calc_category_props(tcomp))
    with open(out_dir / f"{s1}~{s2}.json", 'w') as f:
        json.dump(tcomp, f, sort_keys=False)