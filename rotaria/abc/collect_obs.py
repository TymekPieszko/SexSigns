from pathlib import Path
import json, sys
from rotspec_utils.calc import calc_category_props, rm_silwood, rs_silwood

##############
### Run as:
# python collect_obs.py Rmag silwood 10 15
# python collect_obs.py Rsoc silwood 10 15
# python collect_obs.py Rmag total 10 15
# python collect_obs.py Rsoc total 10 15

##############
### Input:
sp = sys.argv[1]
mode = sys.argv[2]
min_info = int(sys.argv[3])
min_dbl = int(sys.argv[4])
with open(f"../arg_based/kappa/{sp}.json") as f:
    kappa = json.load(f)
with open(f"../arg_based/delta/{sp}/d~{min_dbl}.json") as f:
    delta = json.load(f)
tcomp_dir = Path(f"../arg_based/tcomp/{sp}/i~{min_info}")
silwood = {"Rmag": rm_silwood, "Rsoc": rs_silwood}[sp]

##############
### Output:
out_file = Path(f"./observed/{sp}_{mode}_stats.tsv")
out_file.parent.mkdir(exist_ok=True)

##############
### Main:
# One row per pair
with open(out_file, "w") as f:
    f.write("ind1\tind2\tkappa\tCL\tAR\tSX\tdelta\n")
    for pair in kappa:
        if mode == "silwood" and not set(pair.split("~")) <= set(silwood):
            continue
        with open(tcomp_dir / f"{pair}.json") as g:
            props = calc_category_props(json.load(g))
        f.write("\t".join(map(str, [*pair.split("~"), kappa[pair], props["CL"], props["AR"], props["SX"], delta[pair]])) + "\n")
