from sticcs import sticcs
import numpy as np

def run_sticcs(G, positions, seq_len, forced_ploidy, second_chances=False, silent=False):
    # (TP) From https://github.com/simonhmartin/testing-sticcs-twisst2/blob/main/sim/sim_sticcs.py#L107
    n = G.shape[1]
    # (SM) Check that you have a multiple of forced ploidy haplotypes to work with!
    assert n % forced_ploidy == 0
    ploidies = np.array([forced_ploidy]*int(n/forced_ploidy))
    der_counts, der_positions = sticcs.get_dac_from_haploid_matrix(G, positions, ploidies)
    # print("Derived counts, derived positions:")
    # print(der_counts, der_positions)
    
    # (TP) Example input to the test the above:
    # G = np.array([[0,0,1,1,0,0],[2,2,1,1,3,3],[0,0,1,1,2,2]])
    # get_dac_from_haploid_matrix does something funny with multiallelic sites;
    # e.g., the one above is split into to records in der_counts, 
    # only for both to be removed as containing mising data.  

    # (TP) It could be just:
    # is_bi = np.isin(G, [0, 1]).all(axis=1)
    # G, pos = G[is_bi], pos[is_bi]
    # G = G.reshape(G.shape[0], 2, 2) # Diploidise
    # der_counts = G.sum(axis=2)
    
    # (SM) Currently doesn't accept missing data, so remove any sites with missing data
    no_missing = np.where(der_counts.min(axis=1) >= 0)
    
    der_counts = der_counts[no_missing]
    der_positions = der_positions[no_missing]
    
    # (SM) sticcs has three steps - find patterns, make clusters, infer trees.
    # These could theoretically be combined into a single function within sticcs, but currently I prefer the transparency.
    patterns, matches, n_matches = sticcs.get_patterns_and_matches(der_counts)
    # print("Patterns, matches:")
    # print(patterns, matches)
    
    clusters = sticcs.get_clusters(patterns, matches, der_positions, ploidies=ploidies,
                                    second_chances=second_chances,
                                    seq_start=1, seq_len=seq_len, silent=silent)
    # print("Clusters:")
    # print(clusters)
    
    ts_sticcs = sticcs.infer_ts(patterns, ploidies, clusters, silent=silent)
    
    #simplify is needed for kc_distance to work
    return ts_sticcs #.simplify(reduce_to_site_topology=True)