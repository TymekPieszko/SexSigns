import gzip
from itertools import groupby
import numpy as np
import allel

rm_silwood = ['MAG1', 'MAG2', 'MAG3', 'RMAG402', 'RMAG403', 'RMAG404', 'RMAG407', 'RMAG707'] # Silwood
rm_other = ['RM9', 'RM15'] # Both Italian
rs_silwood = ['AK11', 'AK15', 'AK27', 'RSOC401', 'RSOC405', 'RSOC406', 'RSOC507', 'RSOC702'] # Silwood
rs_other = ['AK16', 'RS1'] # AK16 from Chobham Stream (UK), RS1 - Italian

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

###############
### Rotaria ###
###############

def polarise_genos(genos, alt, anc):
    # Polarise genotypes;
    # returns POLARISED SITES ONLY.
    keep = np.isin(anc, list("ATGC"))
    genos, alt, anc = genos[keep], alt[keep], anc[keep]
    flip = (anc == alt)[:, None, None]
    genos_out = allel.GenotypeArray(np.where(flip, 1 - np.asarray(genos), genos))
    genos_out.is_phased = genos.is_phased
    return genos_out, keep

def get_info_sites(genos, include_sgl=False):
    # Takes 2-ind scikit-allel-style genos;
    # all sites MUST BE polarised;
    # returns INFORMATIVE SITES only.
    # Examples of info sites:
    # 0|1 + 1|0
    # 1/1 + 0|1
    # 0/0 + 1/1
    
    # REMOVED THE GUARD
    # Guard against empty genos;
    # can happen after polarise_genos().
    # if genos.shape[0] == 0:
    #     return genos, np.zeros(0, dtype=bool)

    ac = genos.count_alleles(max_allele=1)
    # Only biallelic with no missingness
    mask_1 = ac.is_biallelic_01() 
    mask_1 = mask_1 & genos.is_called().all(axis=1)
    # Drop sites with an unphased het in either sample
    mask_2 = ~(genos.is_het() & ~genos.is_phased).any(axis=1)
    # Only doubletons and tripletons (singletons too if include_sgl)
    mask_3 = ~ac.is_singleton(allele=1) | include_sgl
    keep = mask_1 & mask_2 & mask_3

    return genos[keep], keep

def get_sgls(genos):
    # Takes 2-ind scikit-allel-style genos;
    # all sites MUST BE polarised;
    # returns SINGLETONS only.
    # Examples of phased singletons:
    # 0|1 + 0|0
    # 0|0 + 1|0
    
    # Guard against empty genos;
    # can happen after polarise_genos().
    if genos.shape[0] == 0:
        return genos, np.zeros(0, dtype=bool)

    ac = genos.count_alleles(max_allele=1)
    # Exactly one derived copy across the 4 haplotypes
    mask_1 = ac.is_singleton(allele=1)
    # Only biallelic with no missingness
    mask_2 = ac.is_biallelic_01()
    mask_2 = mask_2 & genos.is_called().all(axis=1)
    # Drop sites with an unphased het in either sample;
    # a singleton IS a het, so this is what makes it PHASED.
    mask_3 = ~(genos.is_het() & ~genos.is_phased).any(axis=1)
    keep = mask_1 & mask_2 & mask_3

    return genos[keep], keep

def get_phased_blocks(pos, contig_length, ps, is_het):
    is_phased = ps > 0
    
    # Guard against no phased sites
    if sum(is_phased) == 0:
        raise  ValueError("No phased sites!")     

    # Unphased heterozygous SNPs break runs:
    # count them cumulatively, so phased SNPs separated by one get different counts
    is_unphased_het = is_het & ~is_phased
    n_breaks = np.cumsum(is_unphased_het)

    # Contiguous runs of phased SNPs with the same ps tag and break count
    idx = np.flatnonzero(is_phased)
    first, last, n_phased, k = [], [], [], 0
    for _, grp in groupby(zip(ps[idx], n_breaks[idx])):
        n = sum(1 for _ in grp)
        first.append(idx[k])
        last.append(idx[k + n - 1])
        n_phased.append(n)
        k += n
    first = np.array(first)
    last = np.array(last)
    n_phased = np.array(n_phased)

    # Midpoints to the neighbouring heterozygous SNP outside the run
    # (unphased or with a different ps tag), or contig ends
    het_pos = pos[is_het]
    rank = np.cumsum(is_het) - 1 # Rank of each SNP among heterozygous SNPs
    r_first, r_last = rank[first], rank[last]
    last_het = len(het_pos) - 1
    starts = np.where(r_first == 0, 1,
                      (het_pos[np.maximum(r_first - 1, 0)] + het_pos[r_first]) // 2 + 1)
    ends = np.where(r_last == last_het, contig_length,
                    (het_pos[r_last] + het_pos[np.minimum(r_last + 1, last_het)]) // 2)
    return np.column_stack((starts, ends)), n_phased

def get_phased_segments(blocks_1, blocks_2):
    # Takes two (n, 2) arrays of 1-based INCLUSIVE [start, end] blocks for ONE contig,
    # each sorted and non-overlapping (as returned by get_blocks).
    # Returns their pairwise intersections, same format.
    b1 = np.asarray(blocks_1).reshape(-1, 2)
    b2 = np.asarray(blocks_2).reshape(-1, 2)
    out, i, j = [], 0, 0
    while i < len(b1) and j < len(b2):
        start = max(b1[i, 0], b2[j, 0])
        end = min(b1[i, 1], b2[j, 1])
        if start <= end:
            out.append((start, end))
        # advance whichever block ends first
        if b1[i, 1] < b2[j, 1]:
            i += 1
        else:
            j += 1
    return np.array(out, dtype=int).reshape(-1, 2)

def get_contig_lengths(vcf):
    out = {}
    with gzip.open(vcf, "rt") as f:
        for line in f:
            if line.startswith("##contig=<"):
                d = dict(kv.split("=", 1) for kv in line.split("<", 1)[1].rstrip(">\n").split(",") if "=" in kv)
                out[d["ID"]] = int(d["length"])
            elif line.startswith("#CHROM"):
                break
    return out

def read_bed(bed):
    # Produces a numpy structured array - cool!
    return np.genfromtxt(
        bed, dtype=None, encoding=None,
        names=["contig", "start", "end", "name", "n_info"]
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

def calc_kappa(genos):
    # Takes tskit-style genotypes.
    der_counts = genos.sum(axis=1)
    dbl = np.sum(der_counts == 2)
    trp = np.sum(der_counts == 3)
    total = dbl + trp
    if total == 0:
        return np.nan
    return 1 - ((dbl - trp) / total)

def filter_by_dbl(CL_intervals, genos, pos, min_dbl):
    # Takes tskit-style genotypes; keeps intervals with >= min_dbl doubletons.
    is_dbl = genos.sum(axis=1) == 2
    pos = pos[is_dbl]
    return {tree: [(s, e) for s, e in ivs if np.sum((pos >= s) & (pos < e)) >= min_dbl]
            for tree, ivs in CL_intervals.items()}

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
            # counts = np.array([c_0, c_1, c_2, c_3])
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

def get_biallelic_01_allel(genos, pos):
    # Takes and returns scikit-allel-style genotypes.
    # Retains only 0-1 variant sites.
    ac = genos.count_alleles()
    is_bi = ac.is_biallelic_01() # NOTE: this keeps sites with missingness!!
    return genos[is_bi], pos[is_bi]
