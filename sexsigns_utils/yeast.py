import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import tskit, allel, json
from sexsigns_utils.calc import calc_category_props

def read_vcf(vcf):
    callset = allel.read_vcf(str(vcf), fields=["samples", "variants/POS", "calldata/GT"])
    samples = callset["samples"]
    pos = callset["variants/POS"]
    genos = allel.GenotypeArray(callset['calldata/GT'])
    return genos, pos, samples

def empty_ts(sequence_length):
    tables = tskit.TableCollection(sequence_length)
    for _ in range(4):
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAGCLE,
            time=0.0
        )
    return tables.tree_sequence()

def get_features(features, chrom_id):
    # features = pd.read_csv(features, sep="\t")
    features = features.loc[features["chrom"] == f"chrom_{chrom_id}", ["centr_start", "centr_end", "chrom_length"]]
    centr_start, centr_end, chrom_length = features.iloc[0].astype(int)
    return centr_start, centr_end, chrom_length

def get_targets(clade_dict, clade):
    if not isinstance(clade_dict, dict): 
        with open(clade_dict) as f:
            clade_dict = json.load(f)
    targets = clade_dict[clade]
    return targets

def intersect_intervals(intervals_1, intervals_2):
    out = []
    for s1, e1 in intervals_1:
        for s2, e2 in intervals_2:
            if e1 <= s2 or s1 >= e2:
                continue
            out.append((max(s1, s2), min(e1, e2)))
    return np.array(out, dtype=int).reshape(-1, 2)

def get_callable_regions(callable_df, chrom_id, sample):
    callable_df = callable_df[(callable_df["chrom_id"] == chrom_id) & (callable_df["sample"] == sample)]
    callable_regions = callable_df[["start", "end"]].to_numpy(dtype=int) 
    return callable_regions

def get_callable_mask(callable_regions_per_sample, pos):
    mask = []
    # For each sample...
    for regions in callable_regions_per_sample:
        regions = np.array(regions)
        starts = regions[:,0]
        ends = regions[:,1]
        # Do SNPs fall within any of the intervals?
        inside = (pos[:, None] >= starts) & (pos[:, None] < ends)
        # shape: (n_pos, n_intervals)
        mask.append(np.any(inside, axis=1))
    # Masked sites are outside intervals, hence ~ below
    # Masked = True
    return ~np.column_stack(mask) 

def get_callable_length(callable_regions_per_sample):
    out = []
    for regions in callable_regions_per_sample:
        regions = np.array(regions)
        out.append(np.sum(regions[:,1]-regions[:,0]))
    return out

######################
### TCOMP IN SPACE ###
######################

def calc_stcomp(ts, stcomp):
    # stcomp = {code : [] for code in codes}
    for tree in ts.trees():
        if tree.num_roots!=1:
            stcomp["multi_root"].append(tree.interval)
            continue
        rank = tree.rank()
        ##############
        # Sort by indices
        # n4 = np.sort([n for n in tree.nodes()])[4]
        ##############
        # Sort by time
        nodes = list(tree.nodes())
        np.random.shuffle(nodes)
        nodes.sort(key=tree.time)
        n4 = nodes[4]
        ##############
        if rank == (4,1) and set(tree.children(n4)) == set((0,2)):
            stcomp["b0213"].append(tree.interval)
        elif rank == (4,1) and set(tree.children(n4)) == set((1,3)):
            stcomp["b1302"].append(tree.interval)
        elif rank == (4,2) and set(tree.children(n4)) == set((0,3)):
            stcomp["b0312"].append(tree.interval)
        elif rank == (4,2) and set(tree.children(n4)) == set((1,2)):
            stcomp["b1203"].append(tree.interval)
        elif rank == (3,11):
            stcomp["u0123"].append(tree.interval)
        elif rank == (3,8):
            stcomp["u0132"].append(tree.interval)
        elif rank == (3,3):
            stcomp["u2301"].append(tree.interval)
        elif rank == (3,0):
            stcomp["u2310"].append(tree.interval)
        elif rank == (4,0) and set(tree.children(n4)) == set((0,1)):
            stcomp["b0123"].append(tree.interval)
        elif rank == (4,0) and set(tree.children(n4)) == set((2,3)):
            stcomp["b2301"].append(tree.interval)
        elif rank == (3,10):
            stcomp["u0213"].append(tree.interval)
        elif rank == (3,7):
            stcomp["u0312"].append(tree.interval)
        elif rank == (3,9):
            stcomp["u1203"].append(tree.interval)
        elif rank == (3,6):
            stcomp["u1302"].append(tree.interval)
        elif rank == (3,5):
            stcomp["u0231"].append(tree.interval)
        elif rank == (3,4):
            stcomp["u0321"].append(tree.interval)
        elif rank == (3,2):
            stcomp["u1230"].append(tree.interval)
        elif rank == (3,1):
            stcomp["u1320"].append(tree.interval)
        else:
            stcomp["poly"].append(tree.interval)
    return stcomp

###########################
### YEAST VISUALISATION ###
###########################

color_dict = {
    "b0213" : "#377eb8",  # 0
    "b1302" : "#377eb8",  # 1
    "b0312" : "#377eb8",  # 2
    "b1203" : "#377eb8",  # 3
    "u0123" : "#4daf4a",  # 4
    "u0132" : "#4daf4a",  # 5
    "u2301" : "#4daf4a",  # 6
    "u2310" : "#4daf4a",  # 7
    "b0123" : "#4daf4a",  # 8
    "b2301" : "#4daf4a",  # 9
    "u0213" : "#ff7f00",  # 10
    "u0312" : "#ff7f00",  # 11
    "u1203" : "#ff7f00",  # 12
    "u1302" : "#ff7f00",  # 13
    "u0231" : "#ff7f00",  # 14
    "u0321" : "#ff7f00",  # 15
    "u1230" : "#ff7f00",  # 16
    "u1320" : "#ff7f00",  # 17 
    "poly"  : "white",
    "multi_root" : "white",
    "undefined" : "white" 
}

def plot_category_props(ax, tcomp, x, y, font_size):
    cat_props = calc_category_props(tcomp)

    bars = ax.bar(
        range(len(cat_props)),
        list(cat_props.values()),
        color=["#377eb8", "#4daf4a", "#ff7f00"]
    )

    ax.set_ylim(0, 1)

    if font_size:
        for rect, v in zip(bars, cat_props.values()):
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                rect.get_height() + 0.02,
                f"{v:.3f}",
                ha="center",
                va="bottom",
                fontsize=font_size,
                clip_on=False
            )

def plot_grid(data, life_cycle_dict, samples, figsize_scale, font_size):
    color_dict = {"conventional": "#7F9E94", "preferred sexual": "#D26A43","preferred asexual": "#2B2B2B", "not measured": "azure"}
    n = len(samples)

    fig = plt.figure(figsize=(figsize_scale * n, figsize_scale * n))
    gs = gridspec.GridSpec(n, n, figure=fig, wspace=0.25, hspace=0.25)

    axes = [[None] * n for _ in range(n)]

    for i in range(n):
        # for j in range(i+1,n):
        for j in range(i):
            ax = fig.add_subplot(gs[i, j])
            axes[i][j] = ax

            pair_ij = f"{samples[i]}~{samples[j]}"
            pair_ji = f"{samples[j]}~{samples[i]}"
            tcomp = data.get(pair_ij) or data.get(pair_ji)

            if tcomp is not None:
                plot_category_props(ax, tcomp, samples[j], samples[i], font_size)

            ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    label_fs = 36
    # X axis
    bottom_i = n - 1
    for j in range(n - 1):  # columns 0..n-2
        ax = axes[bottom_i][j]
        if ax is not None:
            sample = samples[j]
            ax.set_xlabel(
                sample,
                fontsize=label_fs,
                rotation=66,
                labelpad=6,
                color=color_dict[life_cycle_dict[sample]]
            )
            ax.xaxis.set_label_position("bottom")
    # Y axis
    left_j = 0
    for i in range(1, n):  # rows 1..n-1
        ax = axes[i][left_j]
        if ax is not None:
            sample = samples[i]
            ax.set_ylabel(
                sample,
                fontsize=label_fs,
                rotation=0,
                ha="right",
                va="center",
                labelpad=6,
                color=color_dict[life_cycle_dict[sample]]
            )

    return fig, axes

def calc_spatial_weights(stcomp_lst, tree_codes):  
    # (1) Collect breakpoints
    breakpoints = set()
    for d in stcomp_lst:
        for c in tree_codes:
            for s, e in d.get(c, []):
                breakpoints.add(s)
                breakpoints.add(e)
    breakpoints = sorted(breakpoints)

    blocks = []
    weights = []

    # (2) Iterate blocks
    total = len(stcomp_lst) # For calculating per-block proportions
    for i in range(len(breakpoints) - 1):
        s_i = breakpoints[i]
        e_i = breakpoints[i+1]
        mid = (s_i + e_i) / 2

        # For each dictionary, determine topology at this segment
        counts = {c: 0 for c in tree_codes}
        for d in stcomp_lst:
            for c in tree_codes:
                for s, e in d.get(c, []):
                    if s <= mid < e:
                        counts[c] += 1
                        break
        
        proportions = {c: counts[c] / total for c in tree_codes}

        blocks.append((s_i, e_i))
        weights.append(proportions)

    return blocks, weights

def plot_spatial_weights(blocks, weights, colors, figsize):
    plt.figure(figsize=figsize)
    
    for (s, e), w in zip(blocks, weights):
        bottom = 0
        for code, prop in w.items():
            if prop == 0:
                continue
            # hatch = "xxx" if code == "poly" else ""
            plt.bar(
                x=s,
                height=prop,
                width=e - s,
                bottom=bottom,
                align='edge',
                color=colors[code],
                # hatch=hatch,
                linewidth=0
            )
            bottom += prop

    plt.ylim(0,1)
    # plt.xlabel("Position")
    # plt.ylabel("CL-AR-SX props")