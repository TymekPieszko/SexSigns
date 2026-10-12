from Bio import AlignIO
from pathlib import Path
import sys, gzip 

### About:
# It is assumed that sequence 1 (the reference) is R. magnacalcarata
# and that sequence 2 (the query) is R. socialis.

### Run as:
# python map_sites.py ../03_maf/R_magnacalcarata_RM9~R_socialis_RS1.maf.gz ../04_polar_vcf

target = sys.argv[1]
maf = sys.argv[2]
out_dir = Path(sys.argv[3])

out = open(str(out_dir / f"{target}_map.tsv"), "w")
if target == "Rmag": # Bases in R. magnacalcarata orientation
    out.write("Rmag_contig\tRmag_pos\tRsoc_contig\tRsoc_pos\tRmag_base\tRsoc_base\n")
if target == "Rsoc": # Bases in R. socialis orientation
    out.write("Rsoc_contig\tRsoc_pos\tRmag_contig\tRmag_pos\tRsoc_base\tRmag_base\n")
with gzip.open(maf, "rt") as fh:
    for aln in AlignIO.parse(fh, "maf"):
        # print(aln)
        seq1, seq2 = aln[0], aln[1] # Rmag and Rsoc, respectively
        ann1, ann2 = seq1.annotations, seq2.annotations
        pos1, pos2 = ann1["start"], ann2["start"]

        ################
        ### Rmag frame
        rmag1, rmag2 = str(seq1.seq), str(seq2.seq)
        flip = ann1["strand"] != ann2["strand"]
        ################
        ### Rsoc frame; only needed when that is the target
        if target == "Rsoc":
            rsoc1 = str(seq1.seq.complement()) if flip else rmag1
            rsoc2 = str(seq2.seq.complement()) if flip else rmag2

        contig1, contig2 = seq1.id, seq2.id 
        length1, length2 = ann1["srcSize"], ann2["srcSize"]
        orient1, orient2 = ann1["strand"] == -1, ann2["strand"] == -1

        for i in range(len(rmag1)):
            # complement() leaves "-" untouched, so the gap test is frame-independent
            a, b = rmag1[i], rmag2[i]
            if a != "-" and b != "-":
                # Forward 1-based on each contig's own strand
                fwd1 = length1 - pos1 if orient1 else pos1 + 1
                fwd2 = length2 - pos2 if orient2 else pos2 + 1
                if target == "Rmag":
                    out.write(f"{contig1}\t{fwd1}\t{contig2}\t{fwd2}\t{a.upper()}\t{b.upper()}\n")
                if target == "Rsoc":
                    out.write(f"{contig2}\t{fwd2}\t{contig1}\t{fwd1}\t{rsoc2[i].upper()}\t{rsoc1[i].upper()}\n")
            if a != "-": pos1 += 1
            if b != "-": pos2 += 1
        

