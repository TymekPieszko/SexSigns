from pathlib import Path
import numpy as np
import json, sys, colorsys
import matplotlib.pyplot as plt
from rotspec_utils.calc import rm_silwood, rs_silwood

##############
### Run as:
# python plot_delta_rot.py Rmag silwood 15
# python plot_delta_rot.py Rsoc silwood 15
# python plot_delta_rot.py Rmag total 15
# python plot_delta_rot.py Rsoc total 15

sp = str(sys.argv[1]) # Rmag or Rsoc
mode = str(sys.argv[2]) # silwood or total
min_dbl = int(sys.argv[3])
in_file = f"../delta/{sp}/d~{min_dbl}.json"
plot_file = Path(f"../plots/delta/{sp}.{mode}.d~{min_dbl}.DELTA.png")
plot_file.parent.mkdir(parents=True, exist_ok=True)
if sp=="Rmag":
    color = (75/255, 41/255, 145/255)
if sp=="Rsoc":
    color = (231/255, 90/255, 155/255)
hue, lum, sat = colorsys.rgb_to_hls(*color)
light = colorsys.hls_to_rgb(hue, lum + (1 - lum) * 0.7, sat)

silwood = {"Rmag": rm_silwood, "Rsoc": rs_silwood}[sp]

with open(in_file) as f:
    data = {k: v for k, v in json.load(f).items() if not np.isnan(v)}
is_silwood = np.array([set(k.split("~")) <= set(silwood) for k in data])
data = np.array(list(data.values()))
if mode == "silwood":
    data, is_silwood = data[is_silwood], is_silwood[is_silwood]
mean = np.mean(data)
print(mean)

plt.figure(figsize=(6.4, 6.9))
plt.axes([1.5 / 6.4, 0.4 / 6.9, 4.5 / 6.4, 6 / 6.9])
offset = np.where(is_silwood, -0.1, 0.1) if mode == "total" else 0.0
x = offset + np.random.normal(loc=0.0, scale=0.03, size=len(data))
for sel, fill in [(~is_silwood, light), (is_silwood, color)]:
    plt.plot(x[sel], data[sel], 'o', markersize=16, markerfacecolor=fill, markeredgecolor="black", markeredgewidth=1.2)
plt.plot([-0.22, 0.22], [mean, mean], color="black", linestyle=(0, (1, 0.5)), linewidth=7)
plt.xlim(-0.4, 0.4)
plt.ylim(0, 1)
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False) 
plt.yticks(fontsize=38)
plt.title(f"{sp} d~{min_dbl}", pad=15)
plt.savefig(plot_file, dpi=150)
