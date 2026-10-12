import sys
from pathlib import Path
from rotspec_utils.calc import rm_silwood, rs_silwood

####################
### Run as:
# python collate_stats.py Rmag silwood
# python collate_stats.py Rsoc silwood
# python collate_stats.py Rmag total
# python collate_stats.py Rsoc total

####################
### Inputs:
sp = str(sys.argv[1])
mode = str(sys.argv[2]) # silwood or total
stats_dir = Path(f"./beds_per_pair/{sp}/stats")

####################
### Outputs:
sample_table = Path(f"./beds_per_pair/{sp}/stats/{sp}_sample_stats.{mode}.tsv")
pair_table = Path(f"./beds_per_pair/{sp}/stats/{sp}_pair_stats.{mode}.tsv")

####################
### Collate; sample and pair stats live in their own subdirs
for table, sub_dir in ((sample_table, stats_dir / "sample"),
                       (pair_table, stats_dir / "pair")):
    paths = sorted(sub_dir.glob("*.stats"))
    if not paths:
        print(f"{sub_dir}: no .stats files, skipping", flush=True)
        continue
    lines = [p.read_text().splitlines() for p in paths]
    rows = [row.split("\t") for header, row in lines]
    if mode == "silwood":
        silwood = {"Rmag": rm_silwood, "Rsoc": rs_silwood}[sp]
        n_labels = 2 if sub_dir.name == "pair" else 1
        rows = [r for r in rows if set(r[:n_labels]) <= set(silwood)]

    # Mean row; label columns (e.g. sample, sample_1) are not numeric,
    # so they get "Mean:" / blanks. Decimals are matched to the source.
    mean_row = []
    for k in range(len(rows[0])):
        vals = [r[k] for r in rows]
        try:
            nums = [float(v) for v in vals]
        except ValueError:
            mean_row.append("Mean:" if k == 0 else "")
            continue
        dec = max(len(v.split(".")[1]) if "." in v else 0 for v in vals)
        mean_row.append(f"{sum(nums) / len(nums):.{dec}f}")

    with open(table, "w") as f:
        f.write(lines[0][0] + "\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
        f.write("\t".join(mean_row) + "\n")
    print(f"{table}: {len(rows)} rows + mean", flush=True)
