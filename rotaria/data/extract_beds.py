from pathlib import Path
import numpy as np
from rotspec_utils.calc import get_contig_lengths, get_phased_blocks, get_phased_segments, get_info_sites, polarise_genos
import sys, allel

####################
### About:

# Note that this takes SNP-only VCFs as input;
# invariant sites have played their role in polarisation.
# SNP-only is crucial, since midpoints are between pos
# - this does not differentiate between invariant sites and SNPs.

# Custom phased block / co-phased segment extraction method 
# allows us to use SNP midpoints as boundries; 
# default WhatsHap blocks start / end at SNP positions.

####################
### Run as:
# python extract_beds.py Rsoc AK15 AK16 3

####################
### Inputs:
sp = str(sys.argv[1])
s1 = str(sys.argv[2])
s2 = str(sys.argv[3])
min_info = int(sys.argv[4])
vcf = Path(f"./vcfs_total/{sp}_polar.SNPs.vcf.gz")

####################
### Outputs:
bed = Path(f"./beds_per_pair/{sp}/{s1}~{s2}.bed")
stats_1 = Path(f"./beds_per_pair/{sp}/stats/sample/{s1}.stats")
stats_2 = Path(f"./beds_per_pair/{sp}/stats/sample/{s2}.stats")
stats_12 = Path(f"./beds_per_pair/{sp}/stats/pair/{s1}~{s2}.stats")
bed.parent.mkdir(parents=True, exist_ok=True)
stats_1.parent.mkdir(parents=True, exist_ok=True)
stats_12.parent.mkdir(parents=True, exist_ok=True)

####################
### Parse VCF
# Note genos and ps need to be reordered later
# as scikit-allel reads in VCF order.
callset = allel.read_vcf(
    str(vcf), samples=[s1, s2],
    fields=["samples", "variants/CHROM", "variants/POS", 
            "calldata/GT", "calldata/PS", "variants/ALT", "variants/AA"]
    )
samples = list(callset["samples"])
contigs = callset["variants/CHROM"]
pos = np.asarray(callset["variants/POS"])
genos = allel.GenotypeArray(callset["calldata/GT"])
ps = np.asarray(callset["calldata/PS"]) # ps tags define phased blocks
alt = np.asarray(callset["variants/ALT"])[:,0] # ALT has 3 slots
anc = callset["variants/AA"]
# Parse contig lengths
contig_length_dict = get_contig_lengths(vcf)

####################
### Reorder genos
i, j = [samples.index(s) for s in (s1, s2)]
genos, ps = genos[:, [i, j]], ps[:, [i, j]]
# And phase masks
genos.is_phased = ps > 0

####################
### Compute a list of (contig, start_idx, end_idx)
# This will be used to extract phased blocks per contig
edges = np.flatnonzero(contigs[1:] != contigs[:-1]) + 1
starts = np.concatenate(([0], edges))
ends = np.concatenate((edges, [len(contigs)]))
contig_info = [(contigs[start], start, end) for start, end in zip(starts, ends)]

####################
### Phased blocks (per sample)
block_info = {}
for contig, start_idx, end_idx in contig_info:
    # Sample 1
    if (ps[start_idx:end_idx, 0] > 0).any():
        block_info[(0, contig)] = get_phased_blocks(pos[start_idx:end_idx], contig_length_dict[contig], ps[start_idx:end_idx, 0], genos[start_idx:end_idx, 0].is_het())
    # Sample 2
    if (ps[start_idx:end_idx, 1] > 0).any():
        block_info[(1, contig)] = get_phased_blocks(pos[start_idx:end_idx], contig_length_dict[contig], ps[start_idx:end_idx, 1], genos[start_idx:end_idx, 1].is_het())

####################
### Co-phased segments (per pair)
# Remember:
# Everything already is for just 2 samples!

seg_rows = []
for contig, start_idx, end_idx in contig_info:
    # Guard against contigs with no phased blocks 
    if (0, contig) not in block_info or (1, contig) not in block_info:
        continue
    blocks_s1 = block_info[(0, contig)][0]
    blocks_s2 = block_info[(1, contig)][0]
    # Guard against XXX (is it needed?)
    if len(blocks_s1) == 0 or len(blocks_s2) == 0:
        continue
    # Compute the segments!!
    segments = get_phased_segments(blocks_s1, blocks_s2)

    # Slice per contig, so every index below is contig-local
    pos_ct = pos[start_idx:end_idx]
    genos_ct = genos[start_idx:end_idx]
    alt_ct = alt[start_idx:end_idx]
    anc_ct = anc[start_idx:end_idx]

    for seg_start, seg_end in segments:
        # Compute the idx of where the segments starts / ends
        seg_start_idx, seg_end_idx = np.searchsorted(pos_ct, [seg_start, seg_end + 1])
        # Get just the genos, alt and anc for the segment
        genos_seg = genos_ct[seg_start_idx:seg_end_idx]
        alt_seg = alt_ct[seg_start_idx:seg_end_idx]
        anc_seg = anc_ct[seg_start_idx:seg_end_idx]
        # Count info sites
        n_info = len(get_info_sites(polarise_genos(genos_seg, alt_seg, anc_seg)[0])[0])
        if n_info >= min_info:
            seg_rows.append((contig, seg_start, seg_end, n_info))

####################
### Write bed; 0-based, half-open.
with open(bed, "w") as f:
    for ct, start, end, n_info in seg_rows:
        f.write("%s\t%d\t%d\t%s_%d_%d\t%d\n"
                 % (ct, start - 1, end, ct, start, end, n_info))

####################
### Write stats
# Reference length = all contigs declared in the VCF header.
ref_length = sum(contig_length_dict.values())

### (A) Pair stats
n_seg = len(seg_rows)
# %_ref and mean length of segments with at least this many info sites
thresholds = (min_info, 10, 25, 50, 100)
lengths = [[end - start + 1 for _, start, end, n in seg_rows if n >= t] for t in thresholds]
pct_ref = [100 * sum(l) / ref_length for l in lengths]
mean_length = [sum(l) / len(l) if l else 0.0 for l in lengths]
with open(stats_12, "w") as f:
    f.write("sample_1\tsample_2\tn_segments\t%_ref\t%_ref_>=10\t%_ref_>=25\t%_ref_>=50\t%_ref_>=100"
            "\tmean_length\tmean_length_>=10\tmean_length_>=25\tmean_length_>=50\tmean_length_>=100\n")
    f.write("%s\t%s\t%d" % (s1, s2, n_seg)
            + "".join("\t%.4f" % v for v in pct_ref)
            + "".join("\t%.1f" % v for v in mean_length) + "\n")

### (B) Sample stats
for idx, s, stats in ((0, s1, stats_1), (1, s2, stats_2)):
    blocks = [block_info[(idx, contig)] for contig, _, _ in contig_info if (idx, contig) in block_info]
    n_blocks = sum(len(b) for b, n in blocks)
    block_bp = sum(int((b[:, 1] - b[:, 0] + 1).sum()) for b, n in blocks)
    block_phased = sum(int(n.sum()) for b, n in blocks)
    with open(stats, "w") as f:
        f.write("sample\t%_ref\tn_blocks\tmean_length\tn_phased\n")
        f.write("%s\t%.4f\t%d\t%.1f\t%.1f\n"
                % (s, 100 * block_bp / ref_length, n_blocks,
                   block_bp / n_blocks if n_blocks else 0.0,
                   block_phased / n_blocks if n_blocks else 0.0))
