import sys
from collections import Counter

map_file = sys.argv[1]
snp_file = sys.argv[2]
ortho_file = sys.argv[3] # Both variant and invariant sites
stats_file = sys.argv[4] # e.g. Rmag_biSNPs_polarStats.tsv
out_file = sys.argv[5] # e.g. Rmag_biSNPs_AA.tsv; bgzip + tabix it for bcftools annotate

COMP = str.maketrans("ACGT", "TGCA")

### (1) Dict of focal alleles.
# snp_file can contain biallelic SNPs only;
# others don't matter for our purposes.
alleles = {}
for line in open(snp_file):
    contig, pos, ref, alt = line.rstrip("\n").split("\t")[:4]
    alleles[(contig, int(pos))] = (ref.upper(), alt.upper())

### (2) Dict of orthologous coordinates and bases
sites = {}
with open(map_file) as fh:
    next(fh)
    for line in fh:
        contig, pos, ortho_contig, ortho_pos, _, ortho_base = line.rstrip("\n").split("\t")
        key = (contig, int(pos))
        if key in alleles:
            sites[key] = (ortho_contig, int(ortho_pos), ortho_base.upper())

### (3) Dict of orthologous status, for the coordinates we need only.
need = {(ortho_contig, ortho_pos) for ortho_contig, ortho_pos, _ in sites.values()}
ortho_status = {}
for line in open(ortho_file):
    contig, pos, ref, alt = line.rstrip("\n").split("\t")[:4]
    key = (contig, int(pos))
    if key in need:
        # alt == "." is invariant; the alleles are kept to tell a SHARED
        # polymorphism from an independent one at the same position.
        ortho_status[key] = (ref.upper(), alt.upper())

### (4) Collect stats & write ancestral states (AA) to file
n = Counter()
out = open(out_file, "w")
for key in sorted(alleles):
    if key not in sites:
        n["unaligned"] += 1; continue
    ortho_contig, ortho_pos, ortho_base = sites[key]
    if ortho_base not in alleles[key]:
        n["third allele"] += 1; continue
    record = ortho_status.get((ortho_contig, ortho_pos))
    if record is None: # Absent from the comparator: nothing to polarise against
        n["no comparator record"] += 1; continue
    ortho_ref, ortho_alt = record
    if ortho_alt != ".":
        if ortho_base == ortho_ref:
            flip = False
        elif ortho_base == ortho_ref.translate(COMP):
            flip = True
        else:
            n["comparator REF mismatch"] += 1; continue
        ortho_alleles = {ortho_base} | {a.translate(COMP) if flip else a for a in ortho_alt.split(",")}
        if ortho_alleles == set(alleles[key]):
            n["shared polymorphism"] += 1
        else:
            n["polymorphic, other alleles"] += 1
        continue
    n["AA written"] += 1
    out.write(f"{key[0]}\t{key[1]}\t{ortho_base}\n")
out.close()

### (5) Stats.
total = sum(n.values())
with open(stats_file, "w") as fh:
    fh.write("category\tn\tpct_of_biallelic_SNPs\n")
    for k, v in n.most_common():
        fh.write(f"{k}\t{v}\t{100 * v / total:.3f}\n")
        print(f"{k:28s} {v:>10,}  {100 * v / total:5.1f}%", file=sys.stderr)
    fh.write(f"total_biallelic_SNPs\t{total}\t100.000\n")
print(f"{'total biallelic SNPs':28s} {total:>10,}", file=sys.stderr)
