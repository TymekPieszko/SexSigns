import allel
import numpy as np
import pandas as pd
import itertools
import tqdm
from collections import defaultdict

########################
### ROTARIA ANALYSIS ###
########################
color_rm = (75/255, 41/255, 145/255)
color_rs = (231/255, 90/255, 155/255)
rm_SIL_samples = ['MAG1', 'MAG2', 'MAG3', 'RMAG402', 'RMAG403', 'RMAG404', 'RMAG407', 'RMAG707']
rm_other_samples = ['RM9', 'RM15']
rs_SIL_samples = ['AK11', 'AK15', 'RSOC401', 'RSOC405', 'RSOC406', 'RSOC507', 'RSOC702', 'AK27']
rs_other_samples = ['AK16', 'RS1']

# def get_mask(vcf, ac_min):
#     genos = allel.GenotypeArray(vcf['calldata/GT'])
#     allele_counts = genos.count_alleles()
#     # Biallelic SNP mask
#     bi_mask = np.all(allele_counts!=0, axis=1)
#     # Allele count mask
#     ac_mask = allele_counts[:,1] >= ac_min
#     # Merge masks 
#     mask = bi_mask & ac_mask
#     return mask

def iter_scaffolds(genos, pos, scaffolds):
    start_indices = np.unique(scaffolds, return_index=True)[1]
    scaffold_indices = np.split(np.arange(len(scaffolds)), start_indices[1:])
    for indices in scaffold_indices:
        yield genos[indices], pos[indices]

def iter_ind_pairs(genos, ind_names):
    ind_num = len(ind_names)
    ind_pairs = itertools.product(range(ind_num), repeat=2)
    for i, j in ind_pairs:
        yield (i, j), genos[:,[i,j]]

def calc_heterozygosity_vcf(genos, pos, ind_names, scaffolds):
    # Iterate over i individuals, each with k scaffolds.
    total = defaultdict(dict)
    for i in tqdm.tqdm(range(len(ind_names))):
        genos_i = genos[:,i]
        het_per_individual = []
        for genos_ik, pos_k in iter_scaffolds(genos_i, pos, scaffolds):
            het_num, intervals = allel.windowed_statistic(pos_k, genos_ik, statistic=lambda x : x.is_het(allele=1).sum(), start=0, stop=pos_k[-1], size=5001)[0:2]
            # No snps in an interval --> nan
            het_num = np.nan_to_num(het_num)
            het_per_scaffold = het_num / (intervals[:,1]-intervals[:,0])
            het_per_individual.extend(het_per_scaffold[:-1])
        total[ind_names[i]] = het_per_individual  
    return total

def calc_between_diversity_vcf(genos, pos, scaffolds):
    all_idx = itertools.combinations(np.arange(genos.n_samples*2), r=2)
    # between_idx = [x for x in all_idx if not (x[1] - x[0] == 1 and x[1] % 2 == 1)]
    between_idx = [(i, j) for i in range(0, genos.n_samples*2, 2) for j in (i+1,)]
    total = []
    for genos_k, pos_k in tqdm.tqdm(iter_scaffolds(genos, pos, scaffolds)):
        genos_k = allel.HaplotypeArray(genos_k.reshape(genos_k.shape[0], -1)) # Or rather, 'pseudo haplotypes' or tskit-style genotypes
        genos_k = allel.GenotypeArray([[[row[i], row[j]] for i,j in between_idx] for row in genos_k])
        between_diffs, intervals = allel.windowed_statistic(pos_k, genos_k, statistic=lambda x : x.is_het(allele=1).mean(), start=0, stop=pos_k[-1], size=5001)[0:2]
        between_diffs = np.nan_to_num(between_diffs)
        between_divs = between_diffs / (intervals[:,1]-intervals[:,0])
        total.extend(between_divs[:-1])
    return np.mean(total)

def calc_pi_vcf(genos, pos, scaffolds):
    total = []
    for genos_k, pos_k in tqdm.tqdm(iter_scaffolds(genos, pos, scaffolds)):
        ac_k = genos_k.count_alleles()
        pis = allel.windowed_diversity(pos_k, ac_k, start=0, stop=pos_k[-1], size=5001)[0]
        total.extend(pis[:-1])
    return np.mean(total)



######################
### UNPHASED STATS ###
######################

def get_biallelic(ts):
    genos = allel.HaplotypeArray(ts.genotype_matrix()).to_genotypes(ploidy=2)
    ac = genos.count_alleles()
    bi_mask = ac.is_biallelic_01()
    pos = ts.sites_position
    genos, pos = genos[bi_mask], pos[bi_mask]
    return genos, pos

def calc_delta_branch_based(ts, span_norm):
    """
    Calculates average branch-based delta for individual trees.
    """
    deltas = []
    spans = []
    for tree in ts.trees():
        rank = tree.rank()
        if rank in [(4,1), (3,10), (3,6), (3,5), (3,1)]:
            b1 = tree.time(tree.parent(0)) + tree.time(tree.parent(2))
            b2 = tree.time(tree.parent(1)) + tree.time(tree.parent(3))
            delta = abs(b1-b2) / (b1+b2)
            deltas.append(delta)
            span = tree.span
            spans.append(span)
        elif rank in [(4,2), (3,9), (3,7), (3,4), (3,2)]:
            b1 = tree.time(tree.parent(0)) + tree.time(tree.parent(3))
            b2 = tree.time(tree.parent(1)) + tree.time(tree.parent(2))
            delta = abs(b1-b2) / (b1+b2)
            deltas.append(delta)
            span = tree.span
            spans.append(span)
    if span_norm:
        delta = np.average(deltas, weights=spans)
    else:
        delta = np.mean(deltas)
    span = np.sum(spans)
    return delta, span


def get_tree_intervals(ts, merged):
    """
    Get intervals of incongruent trees; optionally merge intervals of adjacent incongruent trees.
    """
    target_ranks = {(4,1), (3,10), (3,6), (3,5), (3,1), (4,2), (3,9), (3,7), (3,4), (3,2)}
    rank_mask = [tree.rank() in target_ranks for tree in ts.trees()]
    breaks = ts.breakpoints(as_array=True)
    if merged:
        target_switches = np.concatenate(([0], np.where(np.diff(rank_mask))[0]+1, [len(rank_mask)]))
        intervals = [[float(breaks[s]), float(breaks[e])] for s, e in np.column_stack((target_switches[:-1], target_switches[1:]))]
        intervals = intervals[0::2] if rank_mask[0] else intervals[1::2]
        return intervals
    else: 
        intervals = np.array([[float(s), float(e)] for s, e in np.column_stack((breaks[:-1], breaks[1:]))])
        intervals = intervals[rank_mask]
        return intervals

def calc_delta(genos, pos, intervals, mode):
    """
    Calculates average mutation-based delta for given intervals.
    """
    ac = genos.count_alleles()
    if ac.n_alleles < 2:
        return np.nan   
    if mode == "1s_only":
        mask = ac.is_singleton(allele=1)
    elif mode == "3s_explicit" or mode == "3s_as_1s":
        mask = ac.is_singleton(allele=1) | ac.is_singleton(allele=0) # Mask for singletons and tripletons
    genos, pos = genos[mask], pos[mask]
    
    deltas = []
    for start, end in intervals:
        mask_i = (pos>=start) & (pos<end) # Interval mask
        ## Skip empty intervals
        if not np.any(mask_i):
            continue
        genos_1, genos_2 = genos[mask_i,0], genos[mask_i,1]
        if mode == "1s_only":
            c1, c2 = sum(genos_1.is_het()), sum(genos_2.is_het())
        elif mode == "3s_explicit":
            c1, c2 = np.sum(genos_1), np.sum(genos_2)
        elif mode == "3s_as_1s":
            c1 = sum(genos_1.is_het()) + sum(genos_1.is_hom(allele=1))
            c2 = sum(genos_2.is_het()) + sum(genos_2.is_hom(allele=1))
        delta = abs(c1-c2) / (c1 + c2)
        deltas.append(delta)
    delta = np.nanmean(deltas)
    return delta

###################
### Version where the mean is normalised by counts / spans
###################
# def calc_delta(genos, pos, intervals, mode):
#     """
#     Calculates average mutation-based delta for given intervals.
#     """
#     ac = genos.count_alleles()
#     if mode == "1s_only":
#         mask = ac.is_singleton(allele=1)
#     elif mode == "3s_explicit" or mode == "3s_as_1s":
#         mask = ac.is_singleton(allele=1) | ac.is_singleton(allele=0) # Mask for singletons and tripletons
#     genos, pos = genos[mask], pos[mask]
    
#     deltas = []
#     counts = []
#     for start, end in intervals:
#         mask_i = (pos>=start) & (pos<end) # Interval mask
#         ## Skip empty intervals
#         if not np.any(mask_i):
#             continue
#         genos_1, genos_2 = genos[mask_i,0], genos[mask_i,1]
#         if mode == "1s_only":
#             c1, c2 = sum(genos_1.is_het()), sum(genos_2.is_het())
#         elif mode == "3s_explicit":
#             c1, c2 = np.sum(genos_1), np.sum(genos_2)
#         elif mode == "3s_as_1s":
#             c1 = sum(genos_1.is_het()) + sum(genos_1.is_hom(allele=1))
#             c2 = sum(genos_2.is_het()) + sum(genos_2.is_hom(allele=1))
#         delta = abs(c1-c2) / (c1 + c2)
#         deltas.append(delta)
#         count = (c1 + c2)
#         # count = end - start
#         counts.append(count)
#     deltas = np.array(deltas)
#     counts = np.array(counts)
#     mask = ~np.isnan(deltas)
#     deltas = deltas[mask]
#     counts = counts[mask]
#     if counts.sum() == 0:
#         return np.nan
#     delta = np.average(deltas, weights=counts)
#     return delta



def get_double_het_runs(genos, pos, seq_length):
    """
    Get intervals corresponding to runs of double-heterozygous sites.
    """
    ac = genos.count_alleles()
    if ac.n_alleles < 2:
        return [[]]
    is_doubleton = ac.is_doubleton(allele=1)
    if len(is_doubleton) == 0:
        return [[]]
    pos_2 = pos[is_doubleton]
    is_2het_doubleton = np.all(genos[is_doubleton].is_het(), axis=1)
    if len(is_2het_doubleton) == 0:
        return [[]]
    switches = np.where(np.diff(is_2het_doubleton))[0]
    breaks = np.concatenate(([0], ((pos_2[switches]+pos_2[switches+1])/2), [seq_length]))
    intervals = [[float(s), float(e)] for s, e in np.column_stack((breaks[:-1], breaks[1:]))]
    intervals = intervals[0::2] if is_2het_doubleton[0] else intervals[1::2]
    if intervals == []:
        return [[]]
    return intervals

def calc_doubleton_props(genos, ind_names, df=False):
    """
    Returns a 2-D numpy array (or dataframe) with the mean proportion of congruent (homozygous) doubletons for each pair of individuals.
    """
    ind_num = len(ind_names)
    ind_pairs = itertools.product(range(ind_num), repeat=2)
    out = np.zeros(shape=(ind_num, ind_num))
    for i, j in ind_pairs:
        genos_i = genos[:, i]
        genos_j = genos[:, j]
        double_het = sum(genos_i.is_het(1) & genos_j.is_het(1))
        diff_hom = sum(
            genos_i.is_hom(0) & genos_j.is_hom(1)
            | genos_i.is_hom(1) & genos_j.is_hom(0)
        )
        p = diff_hom / (diff_hom + double_het)
        out[i, j] = p
    if df:
        out = pd.DataFrame(out, index=ind_names, columns=ind_names)
    return out

def calc_tripleton_doubleton_ratios(genos, ind_num):
    """
    Returns a 2-D numpy array (or dataframe) with the mean tripleton-to-doubleton ratio for each pair of individuals.
    """
    ind_pairs = itertools.product(range(ind_num), repeat=2)
    out = np.zeros(shape=(ind_num, ind_num))
    for i, j in tqdm.tqdm(ind_pairs):
        genos_ij = genos[:, [i, j]]
        ac = genos_ij.count_alleles()
        t_num = sum(ac.is_singleton(allele=0)) # Tripleton number
        d_num = sum(ac.is_doubleton(allele=1)) # Doubleton number
        # doubletons = sum(genos_ij[:,0].is_het(1) & genos_ij[:,1].is_het(1))
        dt_num = d_num+t_num
        if dt_num == 0:
            ratio = np.nan
        else:
            ratio = 1 - ((d_num-t_num)/(d_num+t_num))
        out[i, j] = ratio
    return out