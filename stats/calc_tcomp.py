from pathlib import Path
import tskit, json, sys
from tqdm import tqdm
from sexsigns_utils.calc import codes, calc_tcomp_vTracked

#############
### ABOUT ###
#############
# Model: GC or CO
# Mode: sim, sticcs or singer
# Bias: e.g., 1.0 or 0.1 

### Run as:
# python calc_tcomp.py GC sim 1.0
# python calc_tcomp.py GC sim 0.05
# python calc_tcomp.py GC sticcs 0.1
# python calc_tcomp.py GC singer 0.1

sim_model = sys.argv[1]
inf_method = sys.argv[2]
if inf_method == "sim":
    target_dir = "1.0.sub_ts"
elif inf_method == "sticcs":
    target_dir = "2.0.sticcs"
elif inf_method == "singer":
    target_dir = "2.1.singer"
else:
    print("Error: unknown inf_method.")
    exit()
bias = float(sys.argv[3])

###################
### Input:
target = Path(f"../sim_pipeline_2/sim_output/{sim_model}/{target_dir}/")
ordered_window_dict = f"./ordered_windows/ordered_windows_{sim_model}.sub_ts.json"
with open(ordered_window_dict, "r") as f:
    ordered_window_dict = json.load(f)
###################
### Output:
out_file = Path(f"./tcomp/{sim_model}_{inf_method}_b~{bias}.json")
out_file.parent.mkdir(exist_ok=True)

total = {}
for dir in target.glob("SEX~*/REC~*/MUT~*"):
    sex = float(dir.parts[-3].split("~")[1])
    rec = float(dir.parts[-2].split("~")[1])
    key = f"{sex}_{rec}"
    print(key)

    # Tree composition will be accumulated across replicates
    tcomp = {code: 0.0 for code in codes}

    ### CASE 1 - ts per simulation! ###
    if inf_method in {"sim", "sticcs"}:

        for ts_file in tqdm(list(dir.glob("*.trees"))):

            # Get ordered_windows
            i = ts_file.stem
            ordered_windows = ordered_window_dict[key][i]
            n_keep = max(1, round(len(ordered_windows) * bias))
            windows = sorted(ordered_windows[:n_keep], key=lambda x: x[0]) # Coordinate sorting required for ts.keep_intervals()
            
            # Calc tcomp
            ts = tskit.load(ts_file)
            ts = ts.keep_intervals(windows).trim()
            tcomp = calc_tcomp_vTracked(ts, tcomp)

    ### CASE 2 - dir per simulation; MCMC samples! ###
    elif inf_method == "singer":

        for ts_dir in tqdm(list(dir.glob("*"))):
            if not ts_dir.is_dir():
                continue
            
            # Get ordered_windows
            i = ts_dir.name
            ordered_windows = ordered_window_dict[key][i]
            n_keep = max(1, round(len(ordered_windows) * bias))
            windows = sorted(ordered_windows[:n_keep], key=lambda x: x[0]) # Coordinate sorting required for ts.keep_intervals()

            # Calc tcomp
            for ts_file in ts_dir.glob("*.trees"):
                # 5 MCMC samples per replicate 
                if int(ts_file.stem.split("_")[1]) % 10 != 0:
                    continue
                ts = tskit.load(ts_file)
                ts = ts.keep_intervals(windows).trim()
                tcomp = calc_tcomp_vTracked(ts, tcomp)

    total[key] = tcomp

with open(out_file, "w") as f:
    json.dump(total, f)