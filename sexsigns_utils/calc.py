import numpy as np
from itertools import groupby
from scipy.stats import norm
import scipy.spatial.distance as ssd
import tskit, allel

############################
### GENERAL / STATISTICS ###
############################

def read_inds(inds_file):
    # Reads an individual file (e.g., inds.txt),
    # where rows are replicates and columns are target individuals.
    # Returns individual pairs per replicate.
    reps = []
    with open(inds_file, "r") as f:
        for line in f:
            inds = list(map(int, line.split("\t")))
            ind_pairs = list(zip(inds[0::2], inds[1::2]))
            reps.append(ind_pairs)
    return reps

def subsample_ts(ts, n):
    # Subsamples a ts to n random individuals.
    inds = np.random.choice(ts.num_individuals, size=n, replace=False)
    nodes = np.concatenate([ts.individual(ind).nodes for ind in inds])
    return ts.simplify(nodes)

def subsample_to_inds(ts, ind_1, ind_2):
    # Subsamples a ts to two given individuals.
    nodes_1 = ts.individual(ind_1).nodes
    nodes_2 = ts.individual(ind_2).nodes
    nodes = list(nodes_1) + list(nodes_2)
    return ts.simplify(samples=nodes)

def get_biallelic_01(genos, pos):
    # Takes and returns tskit-style genotypes.
    # Retains only 0-1 variant sites.
    is_bi = (genos.min(axis=1) == 0) & (genos.max(axis=1) == 1)
    return genos[is_bi], pos[is_bi]

def get_biallelic_01_allel(genos, pos):
    # Takes and returns scikit-allel-style genotypes.
    # Retains only 0-1 variant sites.
    ac = genos.count_alleles()
    is_bi = ac.is_biallelic_01() # NOTE: this keeps sites with missingness!!
    return genos[is_bi], pos[is_bi]

def to_allel(genos):
    # Converts to scikit allel-style genotypes.
    genos = allel.HaplotypeArray(genos).to_genotypes(ploidy=2)
    return genos

def to_tskit(genos):
    # Converts to tskit-style genotypes.
    genos = allel.HaplotypeArray(genos.reshape(genos.shape[0], -1))
    return np.array(genos)

def calc_hi(genos, callable_length):
    # Takes scikit allel-style genotypes.
    # Callable_length can be an integer / float
    # or a 1-D array of length = len(samples).
    
    # Het sites per ind
    if genos.mask is not None:
        het = (genos[:, :, 0] != genos[:, :, 1]) & ~genos.mask
    else:
        het = (genos[:, :, 0] != genos[:, :, 1])
    # Num of hets per ind / callable length per ind
    hi_per_ind = het.sum(axis=0) / callable_length
    hi = np.nanmean(hi_per_ind)
    return hi

########################
### Some alternative ways to calculate heterozygosity
# def calc_hi(ts):
#     hi = ts.diversity(
#         sample_sets=[i.nodes for i in ts.individuals()],
#         mode="site"
#     )
#     hi = np.mean(hi)
#     return hi
# def calc_hi(genos, seq_len):
#     # genos = to_allel(genos)
#     # Het sites per ind
#     het = genos[:, :, 0] != genos[:, :, 1]
#     # Divide by number of sites
#     hi_per_ind = het.sum(axis=0) / seq_len
#     hi = np.mean(hi_per_ind)
#     return hi
########################

def calc_fis(genos):
    # Takes scikit allel-style genotypes.
    fis = np.nanmean(allel.inbreeding_coefficient(genos))
    return fis

def calc_ld(genos, pos, k, d_max):
    # Takes scikit allel-style genotypes.
    n_var = len(genos)
    k = min(k, n_var)

    r2_total = []
    for _ in range(k):
        i = np.random.randint(0, n_var)
        g_i = genos[i]
        p_i = pos[i]
        d = np.abs(pos - p_i)
        valid = (d > 0) & (d < d_max)
        if not np.any(valid):
            continue

        n_i = g_i.to_n_alt()[None, :]
        n_other = genos[valid].to_n_alt()

        r = allel.rogers_huff_r_between(n_i, n_other)
        r2_total.extend(np.ravel(r) ** 2)

    if len(r2_total) == 0:
        return np.nan

    return float(np.nanmean(r2_total))

def calc_fold_ld_decay(genos, pos, seq_len, k, d, w):
    # Computes fold LD decay between windows d bp apart.
    # Takes scikit allel-style genotypes.
    # k = number of windowed comparisons
    # d = distance between windows
    # w = window length

    half_w = w / 2
    mask = (pos >= half_w) & (pos < seq_len - half_w)
    genos, pos = genos[mask], pos[mask]

    ### HELPERS ###
    def get_window_genos(p, half_w):
        # Get genos around the focal snp += half_w
        # and in a window d bp away from the focal snp.
        return genos[(pos >= (p - half_w)) & (pos < (p + half_w))]
    def calc_r2_within(genos):
        n_alt = genos.to_n_alt()
        r = allel.rogers_huff_r(n_alt)
        r2 = r ** 2
        return np.nanmean(r2)
    def calc_r2_between(genos_1, genos_2):
        n_alt_1 = genos_1.to_n_alt()
        n_alt_2 = genos_2.to_n_alt()
        r = allel.rogers_huff_r_between(n_alt_1, n_alt_2)
        r2 = r ** 2
        return np.nanmean(r2)

    decay_lst = []
    count = 0
    while count < k:
        p_1 = np.random.choice(pos) # Choose focal SNP
        genos_1 = get_window_genos(p_1, half_w)
        if p_1 + d + half_w < seq_len:
            p_2 = p_1 + d 
        else:
            p_2 = p_1 - d
        genos_2 = get_window_genos(p_2, half_w)
        # Valid genos
        if len(genos_1) < 2 or len(genos_2) == 0:
            continue

        r2_within = calc_r2_within(genos_1)
        r2_between = calc_r2_between(genos_1, genos_2)
        # Valid r2
        if not np.isfinite(r2_within) or not np.isfinite(r2_between):
            continue
        if r2_between == 0:
            continue

        decay_lst.append(r2_within / r2_between)
        count += 1
    return float(np.nanmean(decay_lst))

def calc_ld_tskit(ts, d1, d2, n_sites=500, seed=None):
    # Takes tskit-style genotypes.

    # Remove sites with !=1 mutation
    ids_to_remove = [site.id for site in ts.sites() if len(site.mutations) != 1]
    ts = ts.delete_sites(ids_to_remove)

    # Get remaining site IDs
    remaining_ids = [site.id for site in ts.sites()]

    # Compute LD
    r2_matrix = tskit.LdCalculator(ts).r2_matrix()
    positions = ts.sites_position
    distances = ssd.squareform(ssd.pdist(positions[:, None]))

    r2 = np.nanmean(r2_matrix[(distances > d1) & (distances < d2)])
    return r2

def calc_r2_windowed(genos, pos, k, window_len, window_num):
    # Takes scikit allel-style genotypes.

    r2_total = []
    dist_total = []

    snps = np.random.choice(len(genos), size=min(k, len(genos)), replace=False)
    for i in snps:
        i = int(i)
        g_i = genos[i]
        p_i = pos[i]

        g_other = allel.GenotypeArray(np.delete(genos, i, axis=0))
        p_other = np.delete(pos, i)

        alt_count_i = [g_i.to_n_alt()]
        alt_counts = g_other.to_n_alt()

        r = np.atleast_1d(
            allel.rogers_huff_r_between(alt_count_i, alt_counts)
        )[0]

        r2 = r ** 2
        dist = np.abs(p_other - p_i)

        mask = (dist != 0) & (dist < window_len * window_num)

        r2_total.extend(r2[mask])
        dist_total.extend(dist[mask])
    r2_total = np.array(r2_total)
    dist_total = np.array(dist_total)

    sort_idx = np.argsort(dist_total)
    r2_total, dist_total = r2_total[sort_idx], dist_total[sort_idx]
    r2 = allel.windowed_statistic(
        dist_total,
        r2_total,
        statistic=lambda x: np.nanmean(x),
        start=0,
        stop=window_len * window_num,
        size=window_len
    )[0]
    return r2

snp_classes = np.array(
    # Classes of doubleton and tripleton SNPs
    [
        [0,0,1,1], # 1
        [1,1,0,0], # 1
        [0,1,0,1], # 2
        [1,0,1,0], # 2
        [0,1,1,0], # 3
        [1,0,0,1], # 3
        [1,1,1,0],
        [1,1,0,1],
        [1,0,1,1],
        [0,1,1,1],
    ]
)

def get_tree_intervals(ts, include_SX):
    # If include_SX == False, each interval contains CL trees only.
    # If include_SX == True, an interval can contain a mix of CL and SX trees
    # sharing the same unrooted topology.
    if include_SX:
        set_1 = {(4,1), (3,10), (3,6), (3,5), (3,1)}
        set_2 = {(4,2), (3,9), (3,7), (3,4), (3,2)}
    else:
        set_1 = {(4,1)}
        set_2 = {(4,2)}

    def classify(tree):
        if tree.num_roots != 1: # Cannot be ranked
            return 0
        rank = tree.rank()
        if rank in set_1:
            return 1
        elif rank in set_2:
            return 2
        else:
            return 0

    trees = (classify(tree) for tree in ts.trees())
    breaks = ts.breakpoints(as_array=True)

    i = 0
    interval_dict = {1: [], 2: []}

    for k, v in groupby(trees):
        increment = sum(1 for _ in v)
        interval = (breaks[i], breaks[i + increment])
        i += increment
        if k in (1, 2):
            interval_dict[k].append(interval)

    return interval_dict

def get_snp_intervals(genos, pos, seq_len, include_trp):
    # If include_trp == False, each interval contains doubletons only.
    # If include_trp == True, an interval can contain a mix of compatible doubletons and tripletons.
    mask = genos.sum(axis=1) == 1
    if include_trp:
        mask |= genos.sum(axis=1) == 3

    genos, pos = genos[~mask], pos[~mask]

    if len(genos) == 0:
        return {1: [], 2: []}

    genos = np.where(
        np.all(genos[:, None, :] == snp_classes[None, :, :], axis=2)
    )[1] // 2

    breaks = np.concatenate(([0.0], (pos[:-1] + pos[1:]) / 2, [seq_len]))

    i = 0
    interval_dict = {1: [], 2: []}

    for k, v in groupby(genos):
        increment = sum(1 for _ in v)
        interval = (breaks[i], breaks[i + increment])
        i += increment
        if k in (1, 2):
            interval_dict[k].append(interval)

    return interval_dict

def calc_delta(genos, pos, interval_dict):
    # Takes tskit-style genotypes.
    is_sgl = genos.sum(axis=1) == 1
    genos, pos = genos[is_sgl], pos[is_sgl]

    # Delta value will be normalised by interval span
    weighted_sum = 0.0
    total_span = 0.0

    for tree, intervals in interval_dict.items():
        for s, e in intervals:

            mask_i = (pos >= s) & (pos < e)
            if sum(mask_i) == 0:
                continue

            genos_i = genos[mask_i]
            c_0, c_1, c_2, c_3 = genos_i.sum(axis=0)

            # To avoid a bias towards near-0 values:
            counts = np.array([c_0, c_1, c_2, c_3])
            # if np.count_nonzero(counts) < 2:
            #     continue 
            
            tree = int(tree)
            if tree == 1:
                w_1 = abs(c_0 - c_2); w_2 = abs(c_1 - c_3)
                b_1 = abs(c_0 - c_1); b_2 = abs(c_0 - c_3)
                b_3 = abs(c_2 - c_1); b_4 = abs(c_2 - c_3)

            if tree == 2:
                w_1 = abs(c_0 - c_3); w_2 = abs(c_1 - c_2)
                b_1 = abs(c_0 - c_1); b_2 = abs(c_0 - c_2)
                b_3 = abs(c_3 - c_1); b_4 = abs(c_3 - c_2)

            total = genos_i.sum()
            if total > 0:
                delta_i = ((b_1 + b_2 + b_3 + b_4) - 2 * (w_1 + w_2)) / total
            else:
                continue

            span = e - s
            weighted_sum += delta_i * span
            total_span += span

    if total_span > 0:
        delta = weighted_sum / total_span
        return delta, total_span
    else:
        return np.nan, 0.0

def count_sgl(genos, pos, interval_dict):
    # Mirroring the processing inside calc_delta()
    is_sgl = genos.sum(axis=1) == 1
    genos, pos = genos[is_sgl], pos[is_sgl]

    sgl_counts = []
    for tree, intervals in interval_dict.items():
        for s, e in intervals:

            mask_i = (pos >= s) & (pos < e)
            if sum(mask_i) == 0:
                continue

            genos_i = genos[mask_i]
            c_0, c_1, c_2, c_3 = genos_i.sum(axis=0)

            counts = np.array([c_0, c_1, c_2, c_3])
            if np.count_nonzero(counts) < 2:
                continue 

            count = np.sum(counts)
            sgl_counts.append(count)

    return np.nanmean(sgl_counts)

def calc_kappa(genos):
    # Takes tskit-style genotypes.
    der_counts = genos.sum(axis=1)
    dbl = np.sum(der_counts == 2)
    trp = np.sum(der_counts == 3)
    total = dbl + trp
    if total == 0:
        return np.nan
    return 1 - ((dbl - trp) / total)

def calc_akaike_weights(obs, dists, n_params):
    # dists should be a dictionary of simulated ditributions
    log_likelihoods = []
    for dist in dists:
        dist = np.array(dist)
        mean = np.nanmean(dist)
        std = np.nanstd(dist)
        likelihood = norm.pdf(obs, loc=mean, scale=std)
        log_likelihoods.append(np.log(likelihood) if likelihood > 0 else -np.inf)

    log_likelihoods = np.array(log_likelihoods)
    aic = 2 * n_params - 2 * log_likelihoods
    delta_aic = aic - np.min(aic)
    rel_likelihoods = np.exp(-0.5 * delta_aic)
    return rel_likelihoods / np.sum(rel_likelihoods)

########################
### TREE COMPOSITION ###
########################

codes = [
    "b0213",  # 0
    "b1302",  # 1
    "b0312",  # 2
    "b1203",  # 3
    "u0123",  # 4
    "u0132",  # 5
    "u2301",  # 6
    "u2310",  # 7
    "b0123",  # 8
    "b2301",  # 9
    "u0213",  # 10
    "u0312",  # 11
    "u1203",  # 12
    "u1302",  # 13
    "u0231",  # 14
    "u0321",  # 15
    "u1230",  # 16
    "u1320",  # 17
]

def calc_tcomp_vRanks(ts, tcomp):
    # tcomp should be a dictionary 
    # e.g., {code : 0.0 for code in codes}
    for tree in ts.trees():
        if tree.num_roots!=1:
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
            tcomp["b0213"] += tree.span
        elif rank == (4,1) and set(tree.children(n4)) == set((1,3)):
            tcomp["b1302"] += tree.span
        elif rank == (4,2) and set(tree.children(n4)) == set((0,3)):
            tcomp["b0312"] += tree.span
        elif rank == (4,2) and set(tree.children(n4)) == set((1,2)):
            tcomp["b1203"] += tree.span
        elif rank == (3,11):
            tcomp["u0123"] += tree.span
        elif rank == (3,8):
            tcomp["u0132"] += tree.span
        elif rank == (3,3):
            tcomp["u2301"] += tree.span
        elif rank == (3,0):
            tcomp["u2310"] += tree.span
        elif rank == (4,0) and set(tree.children(n4)) == set((0,1)):
            tcomp["b0123"] += tree.span
        elif rank == (4,0) and set(tree.children(n4)) == set((2,3)):
            tcomp["b2301"] += tree.span
        elif rank == (3,10):
            tcomp["u0213"] += tree.span
        elif rank == (3,7):
            tcomp["u0312"] += tree.span
        elif rank == (3,9):
            tcomp["u1203"] += tree.span
        elif rank == (3,6):
            tcomp["u1302"] += tree.span
        elif rank == (3,5):
            tcomp["u0231"] += tree.span
        elif rank == (3,4):
            tcomp["u0321"] += tree.span
        elif rank == (3,2):
            tcomp["u1230"] += tree.span
        elif rank == (3,1):
            tcomp["u1320"] += tree.span
    return tcomp

def calc_tcomp_vTracked(ts, tcomp):
    # tcomp should be a dictionary 
    # e.g., {code : 0.0 for code in codes}
    sample_sets = [[0, 1], [2, 3], [0, 2], [1, 3], [0, 3], [1, 2]]
    tracked = [tskit.Tree(ts, tracked_samples=samples) for samples in sample_sets]

    for i in range(ts.num_trees):
        for tree in tracked:
            tree.next()
            assert i == tree.index

        # Ignore trees with polytomies (retain trees with unary nodes)
        if any(
            tree.num_children(u) > 2
            for u in tracked[0].nodes()
            if tracked[0].is_internal(u)
        ):
            continue
        # Ignore trees with multiple roots
        if not tracked[0].has_single_root:
            continue

        nodes = list(tracked[0].nodes(order="timeasc"))
        assert nodes[-1] == tracked[0].root
        node4, node5 = nodes[4:6]

        counts = [
            [tree.num_tracked_samples(node4), tree.num_tracked_samples(node5)]
            for tree in tracked
        ]
        counts_set = [set(c) for c in counts]

        # BALANCED
        if counts_set[0] == counts_set[1] == {0, 2}:  # b0123
            if tracked[0].time(node4) == tracked[0].time(node5):
                # tcomp["b4_eq"] += tracked[0].span
                tcomp["b0123"] += tracked[0].span / 2
                tcomp["b2301"] += tracked[0].span / 2
            else:
                if counts[0][0] == 2:
                    tcomp["b0123"] += tracked[0].span
                if counts[0][0] == 0:
                    tcomp["b2301"] += tracked[0].span
        if counts_set[2] == counts_set[3] == {0, 2}:  # b0213
            if tracked[0].time(node4) == tracked[0].time(node5):
                # tcomp["b5_eq"] += tracked[0].span
                tcomp["b0213"] += tracked[0].span / 2
                tcomp["b1302"] += tracked[0].span / 2
            else:
                if counts[2][0] == 2:
                    tcomp["b0213"] += tracked[0].span
                if counts[2][0] == 0:
                    tcomp["b1302"] += tracked[0].span
        if counts_set[4] == counts_set[5] == {0, 2}:  # b0312
            if tracked[0].time(node4) == tracked[0].time(node5):
                # tcomp["b6_eq"] += tracked[0].span
                tcomp["b0312"] += tracked[0].span / 2
                tcomp["b1203"] += tracked[0].span / 2
            else:
                if counts[4][0] == 2:
                    tcomp["b0312"] += tracked[0].span
                if counts[4][0] == 0:
                    tcomp["b1203"] += tracked[0].span
        # UNBALANCED
        if counts_set[0] == {2}:  # u01XX
            if counts_set[2] == counts_set[5] == {1, 2}:
                tcomp["u0123"] += tracked[0].span
            if counts_set[3] == counts_set[4] == {1, 2}:
                tcomp["u0132"] += tracked[0].span
        if counts_set[1] == {2}:  # u23XX
            if counts_set[2] == counts_set[4] == {1, 2}:
                tcomp["u2301"] += tracked[0].span
            if counts_set[3] == counts_set[5] == {1, 2}:
                tcomp["u2310"] += tracked[0].span
        if counts_set[2] == {2}:  # u02XX
            if counts_set[0] == counts_set[5] == {1, 2}:
                tcomp["u0213"] += tracked[0].span
            if counts_set[1] == counts_set[4] == {1, 2}:
                tcomp["u0231"] += tracked[0].span
        if counts_set[3] == {2}:  # u13XX
            if counts_set[0] == counts_set[4] == {1, 2}:
                tcomp["u1302"] += tracked[0].span
            if counts_set[1] == counts_set[5] == {1, 2}:
                tcomp["u1320"] += tracked[0].span
        if counts_set[4] == {2}:  # u03XX
            if counts_set[0] == counts_set[3] == {1, 2}:
                tcomp["u0312"] += tracked[0].span
            if counts_set[1] == counts_set[2] == {1, 2}:
                tcomp["u0321"] += tracked[0].span
        if counts_set[5] == {2}:  # u12XX
            if counts_set[0] == counts_set[2] == {1, 2}:
                tcomp["u1203"] += tracked[0].span
            if counts_set[1] == counts_set[3] == {1, 2}:
                tcomp["u1230"] += tracked[0].span
    return tcomp

def calc_category_props(tcomp):
    CL = sum(
        v for k, v in tcomp.items() if k in ["b0213", "b1302", "b0312", "b1203", "b5_eq", "b6_eq"]
    )

    AR = sum(
        v for k, v in tcomp.items() if k in ["u0123", "u0132", "u2301", "u2310", "b0123", "b2301", "b4_eq"]
    )
    SX = sum(
        v
        for k, v in tcomp.items()
        if k in ["u0213", "u0231", "u1302", "u1320", "u0312", "u0321", "u1203", "u1230"]
    )
    total = CL + AR + SX
    if total == 0.0:
        return {"CL":0.0, "AR":0.0, "SX":0.0}
    CL = CL / total
    AR = AR / total
    SX = SX / total
    return {"CL":CL, "AR":AR, "SX":SX}