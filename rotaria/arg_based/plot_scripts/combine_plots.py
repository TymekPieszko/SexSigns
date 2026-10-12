from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import sys

##############
### Run as:
# python combine_plots.py silwood 10 15
# python combine_plots.py total 10 15

mode = sys.argv[1] # silwood or total
min_info = int(sys.argv[2])
min_dbl = int(sys.argv[3])
plot_file = Path(f"../plots/Rotaria.{mode}.i~{min_info}.d~{min_dbl}.COMBINED.png")

##############
### Rows: species; columns: kappa, tcomp, delta
grid = [[mpimg.imread(f"../plots/kappa/{sp}.{mode}.KAPPA.png"),
         mpimg.imread(f"../plots/tcomp/{sp}.{mode}.i~{min_info}.TOTAL.png"),
         mpimg.imread(f"../plots/delta/{sp}.{mode}.d~{min_dbl}.DELTA.png")] for sp in ["Rmag", "Rsoc"]]

widths = [img.shape[1] / img.shape[0] for img in grid[0]]
wspace = 0.05
fig, axs = plt.subplots(2, 3, figsize=(4 * (sum(widths) + 2 * wspace * sum(widths) / 3), 8),
                        gridspec_kw={"width_ratios": widths, "wspace": wspace, "hspace": 0})
for row, imgs in zip(axs, grid):
    for ax, img in zip(row, imgs):
        ax.imshow(img)
        ax.axis("off")
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
fig.savefig(plot_file, dpi=300)
