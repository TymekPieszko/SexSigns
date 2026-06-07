import shutil
from pathlib import Path
from functions_smk import run_singer, run_convert_to_tskit
import time

vcf = Path(snakemake.input.vcf)
ts_dir = Path(snakemake.output.ts_dir)
arg_dir = ts_dir / "_arg"
L = str(snakemake.params.L)

# Create output dirs
ts_dir.mkdir(parents=True, exist_ok=True)
arg_dir.mkdir(parents=True, exist_ok=True)


# Params for singer_master
singer_path = "singer_master"
Ne = "1000"
m = "5.0e-7"
n = "100"
thin = "20"
polar = "0.99" # VCO is correctly polarised

# Params for convert_to_ts
convert_to_tskit_path = "convert_to_tskit"
start_index = "50"
end_index = "100"
step = "5"

rep = vcf.stem
vcf_prefix = str(vcf.parent / rep)
arg_prefix = str(arg_dir / rep)
ts_prefix = str(ts_dir / rep)

t_start = time.time() 

run_singer(
    singer_path=singer_path,
    Ne=Ne,
    m=m,
    vcf_prefix=vcf_prefix,
    arg_prefix=arg_prefix,
    start="0",
    end=L,
    n=n,
    thin=thin,
    polar=polar,
)
run_convert_to_tskit(
    convert_to_tskit_path=convert_to_tskit_path,
    arg_prefix=arg_prefix,
    ts_prefix=ts_prefix,
    start=start_index,
    end=end_index,
    step=step,
)

shutil.rmtree(arg_dir)

t_end = time.time()

elapsed = t_end - t_start
time_file = ts_dir / "time.txt"
with open(time_file, "w") as f:
    f.write(str(elapsed))

