from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

##############
### Run as:
# python combine_pca_plots.py

plot_file = Path("pca_combined.png")

##############
### Columns: Silwood-only, all samples
imgs = [mpimg.imread(f"pca_{mode}.png") for mode in ["silwood", "total"]]

widths = [img.shape[1] / img.shape[0] for img in imgs]
wspace = 0.1
fig, axs = plt.subplots(1, 2, figsize=(8 * (sum(widths) + wspace * sum(widths) / 2), 8),
                        gridspec_kw={"width_ratios": widths, "wspace": wspace})
for ax, img in zip(axs, imgs):
    ax.imshow(img)
    ax.axis("off")
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
fig.savefig(plot_file, dpi=300)
