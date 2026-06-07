import tskit
from sticcs import sticcs

def run_sticcs(G, positions, seq_len, forced_ploidy, second_chances, silent):
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

ts_in = snakemake.input[0]
ts_out = snakemake.output[0]
forced_ploidy = int(snakemake.params.forced_ploidy)

def get_sticcs_inputs(ts):
    ts = tskit.load(ts)
    G = ts.genotype_matrix()
    pos = ts.sites_position
    if pos[0] == 0.0:
        G = G[1:]
        pos = pos[1:]
    seq_len = ts.sequence_length
    # Biallelic SNPs only.
    is_bi_01 = [{0,1} == set(g) for g in G]
    G, pos = G[is_bi_01], pos[is_bi_01]
    return G, pos, seq_len

G, pos, seq_len = get_sticcs_inputs(ts_in)
# Run sticcs
ts = run_sticcs(G, pos, seq_len, forced_ploidy=forced_ploidy, second_chances=True, silent=False)
ts.dump(ts_out)