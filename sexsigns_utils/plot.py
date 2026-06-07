import matplotlib.pyplot as plt
import seaborn as sns

# Sex labels
sex_labels = ["0.0", "1.0e-05", "3.162e-05", "1.0e-04", "3.162e-04", "1.0e-03", "3.162e-03", "1.0e-02", "3.162e-02", "0.1", "3.162e-01"]
sex_labels_sparse = ["0.0", "1.0e-05", "1.0e-04", "1.0e-03", "1.0e-02", "0.1"]
sex_labels_in_N = ["0", r"0.01 $N^{-1}$", r"0.03162 $N^{-1}$", r"0.1 $N^{-1}$", r"0.3162 $N^{-1}$", r"$N^{-1}$", r"3.162 $N^{-1}$",  r"10 $N^{-1}$", r"31.62 $N^{-1}$", r"100 $N^{-1}$", r"316.2 $N^{-1}$"]
sex_labels_in_N_sparse = ["0", r"0.01 $N^{-1}$", r"0.1 $N^{-1}$", r"$N^{-1}$", r"10 $N^{-1}$", r"100 $N^{-1}$"]


# Rec labels
rec_labels_GC = ["0.0", "1.0e-09", "3.162e-09", "1.0e-08", "3.162e-08", "1.0e-07", "3.162e-07", "1.0e-06"]
rec_labels_CO = [0.0, 2.0193506995300517e-11, 6.324221574335387e-11, 2.0004670849692646e-10, 6.326010734420079e-10, 2.0023521210661537e-09, 6.344556165159151e-09, 2.020253446031858e-08]
loh_labels = ["0.0", "5.0e-06", "1.58e-05", "5.0e-05", "1.58e-04", "5.0e-04", "1.58e-03", "5.0e-03"]
loh_labels_sparse = ["0.0", "5.0e-06", "5.0e-05", "5.0e-04", "5.0e-03"]
loh_labels_in_N = ["0", r"0.005 $N^{-1}$", r"0.0158 $N^{-1}$", r"0.05 $N^{-1}$", r"0.158 $N^{-1}$", r"0.5 $N^{-1}$", r"1.58 $N^{-1}$", r"5 $N^{-1}$"]
loh_labels_in_N_sparse = ["0", r"0.005 $N^{-1}$", r"0.05 $N^{-1}$", r"0.5 $N^{-1}$", r"5 $N^{-1}$"]


def plot_heatmap(
	data,
    cmap,
    vmin,
    vmax,
    xlabel,
    ylabel,
    tick_labels_in_N,
    title,
    title_font,
    label_font,
    tick_font,
    cbar,
    annot
 ):
    
    if cbar == True:
        figsize = (14, 10)
    else:
        figsize = (12, 10) 
    fig, ax = plt.subplots(figsize=figsize)

    hm = sns.heatmap(
        data,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        annot=annot,
        fmt=".3f",
        ax=ax,
        cbar=cbar,
    )

    plt.xlabel(xlabel, fontsize=label_font, labelpad=-6)
    plt.ylabel(ylabel, fontsize=label_font, labelpad=-6)
    if tick_labels_in_N:
        # xtick_labels = ["0", r"0.01 $N^{-1}$", r"0.1 $N^{-1}$", r"$N^{-1}$", r"10 $N^{-1}$", r"100 $N^{-1}$"]
        # xtick_labels = ["0.0"] + sex_labels_in_N[1:][::2]
        xtick_labels = sex_labels_in_N_sparse
        # ytick_labels = [r"5 $N^{-1}$", r"0.5 $N^{-1}$", r"0.05 $N^{-1}$", r"0.005 $N^{-1}$", "0"]
        ytick_labels = loh_labels_in_N_sparse[::-1]
    else:
        # xtick_labels = ["0.0", "1e-05", "1e-04", "1e-03", "1e-02", "1e-01"]
        # xtick_labels = ["0.0"] + sex_labels[1:][::2]
        xtick_labels = sex_labels_sparse
        # ytick_labels = ["5.0e-03", "5.0e-04", "5.0e-05", "5.0e-06", "0.0"]
        # ytick_labels = loh_labels[1:][::2][::-1] + ["0.0"]
        ytick_labels = loh_labels_sparse[::-1]
    plt.xticks(
        ticks=[0.5, 1.5, 3.5, 5.5, 7.5, 9.5],
        labels=xtick_labels,
        fontsize=tick_font,
        rotation=45,
    )
    plt.yticks(
        ticks=[0.5, 2.5, 4.5, 6.5, 7.5],
        labels=ytick_labels,
        fontsize=tick_font,
        rotation=0,
    )
    if cbar:
        cbar_obj = hm.collections[0].colorbar
        cbar_obj.ax.tick_params(labelsize=tick_font)
    plt.title(title, fontsize=title_font, pad=16)
    plt.tight_layout()
    return fig, ax