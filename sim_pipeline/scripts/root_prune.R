library(ape)
library(phangorn)

fasta_dir <- snakemake@input[[1]]
tree_dir_in <- snakemake@input[[2]]
tree_dir_out <- snakemake@output[[1]]

inds_file <- paste0(fasta_dir, "/inds.txt")
inds <- scan(inds_file, what = numeric(), sep = "\t")
inds <- inds[!is.na(inds)] # Sadly, there are trailing tabs...
# cat(inds)
in_files <- list.files(tree_dir_in, pattern = "*.treefile", full.names = TRUE)
cat(in_files)

dir.create(tree_dir_out)


for (in_file in in_files) {
    # print(in_file)
    tree <- read.tree(file = in_file)
    # print(tree)
    tree <- midpoint(tree)

    # inds <- sample(0:99, size = 2, replace = FALSE)
    samples <- unlist(lapply(inds, function(i) c(paste0("n", i * 2), paste0("n", i * 2 + 1))))
    # print(in_file)
    print(samples)
    tree <- keep.tip(tree, samples)

    out_file <- paste0(tree_dir_out, "/", basename(in_file))
    write.tree(tree, file = out_file)
}
# newick_s <- commandArgs(trailingOnly = TRUE)[1]
# tree <- read.tree(file = newick_s)
# tree <- midpoint(tree)

# inds <- sample(0:99, size = 2, replace = FALSE)
# samples <- unlist(lapply(inds, function(i) c(paste0("n", i * 2), paste0("n", i * 2 + 1))))
# tree <- keep.tip(tree, samples)
# # tree <- chronos(tree, quiet = TRUE)

# cat(write.tree(tree))
