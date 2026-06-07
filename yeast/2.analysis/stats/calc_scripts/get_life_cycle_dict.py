from pathlib import Path
import pandas as pd
import numpy as np
import json

# python get_life_cycle_dict.py

##############
### Input:
clade_dict = f"/data/biol-bdelloids/scro4331/SexSigns_2025/yeast/1.args/sticcs/clades_PH.json"
with open(clade_dict, 'r') as f:
    clade_dict = json.load(f)
df = pd.read_excel(open("Becerra-Rodriguez_et_al_2026.ign.xlsx", "rb"), sheet_name="Sheet 1", header=2) # Supplementary table from the paper

##############
### Output:
out_file = Path("./life_cycle_dict.json")

samples = np.concatenate([x for x in clade_dict.values()])

out_dict = {}
for sample in samples:
    life_cycle = df.loc[df["StandardizedName"] == sample, "Life cycle type"].iloc[0]
    out_dict[sample] = life_cycle

with open(out_file, 'w') as f:
    json.dump(out_dict, f)