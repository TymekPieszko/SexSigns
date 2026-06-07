from pathlib import Path
import pandas as pd
from tqdm import tqdm
import sys, json

# python collect_stats.py GC

###################
### Input:
sim_model = str(sys.argv[1])
hi_data = f"../stats/H_i/{sim_model}.json"
fis_data = f"../stats/F_is/{sim_model}.json"
ld_data = f"../stats/LD_0-1000/{sim_model}.json"
kappa_data = f"../stats/kappa/{sim_model}.json"
delta_bias_data = f"../stats/delta/{sim_model}_dbl_trp_b~0.1.plot.json"
delta_data = f"../stats/delta/{sim_model}_dbl_trp_b~1.0.plot.json"
with open(hi_data, "r") as f:
    hi_data = json.load(f)
with open(fis_data, "r") as f:
    fis_data = json.load(f)
with open(ld_data, "r") as f:
    ld_data = json.load(f)
with open(kappa_data, "r") as f:
    kappa_data = json.load(f)
with open(delta_bias_data, "r") as f:
    delta_bias_data = json.load(f)
with open(delta_data, "r") as f:
    delta_data = json.load(f)

###################
### Output:
out_file = Path(f"./out/{sim_model}_stats.tsv")
out_file.parent.mkdir(exist_ok=True)

###################
### Main:

sex_lst = sorted({float(k.split("_")[0]) for k in hi_data.keys()})
sex_lst = [sex for sex in sex_lst if sex != 1.0]
rec_lst = sorted({float(k.split("_")[1]) for k in hi_data.keys()})
print(sex_lst)
print(rec_lst)

table = []
for sex in tqdm(sex_lst):
    for rec in rec_lst:
        a = hi_data[f"{sex}_{rec}"]
        b = fis_data[f"{sex}_{rec}"]
        c = ld_data[f"{sex}_{rec}"]
        d = kappa_data[f"{sex}_{rec}"]
        e = delta_bias_data[f"{sex}_{rec}"]
        f = delta_data[f"{sex}_{rec}"]
        for idx in range(100):
            table.append([sex, rec, a[idx], b[idx], c[idx], d[idx], e[idx], f[idx]])

cols = ["sex", "rec", "H_i", "F_is", "r2", "kappa", "delta_bias", "delta"]
df = pd.DataFrame(table, columns=cols)
df.to_csv(out_file, sep='\t', index=False)
