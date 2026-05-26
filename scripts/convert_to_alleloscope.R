library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
input_dir  <- args[which(args=="--input-dir")  + 1]
output_dir <- args[which(args=="--output-dir") + 1]

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(paste0(output_dir,"/plots"), showWarnings=FALSE)
dir.create(paste0(output_dir,"/rds"), showWarnings=FALSE)

cat("Loading dataset...\n")
cells   <- read.table(paste0(input_dir,"/cells.tsv"), sep="\t", header=TRUE)
bins    <- read.table(paste0(input_dir,"/bins.tsv"),  sep="\t", header=TRUE)
rc      <- read.table(paste0(input_dir,"/readcounts.tsv"), sep="\t", header=TRUE, row.names=1)
alt_bin <- read.table(paste0(input_dir,"/allele_alt.tsv"), sep="\t", header=TRUE, row.names=1)
ref_bin <- read.table(paste0(input_dir,"/allele_ref.tsv"), sep="\t", header=TRUE, row.names=1)

n_cells    <- nrow(cells)
n_bins     <- nrow(bins)
cell_names <- cells$cell
normal_cells <- cells$cell[cells$clone_assignment == 0]
cat("  Cells:", n_cells, "| Bins:", n_bins, "| Normal:", length(normal_cells), "\n")

# ── Tumor and normal matrices ─────────────────────────────────────────────────
row_names  <- paste0(bins$chrom, "-", bins$start, "-", bins$end)
tumor_mat  <- as.data.frame(t(rc))
rownames(tumor_mat) <- row_names
colnames(tumor_mat) <- cell_names

if (length(normal_cells) > 0) {
    normal_vec <- round(colMeans(as.matrix(rc[normal_cells, , drop=FALSE])))
} else {
    normal_vec <- round(colMeans(as.matrix(rc)))
}
normal_mat <- data.frame(normal=pmax(as.integer(normal_vec), 1L), row.names=row_names)

write.table(tumor_mat,  paste0(output_dir,"/tumor.txt"),  sep="\t", quote=FALSE)
write.table(normal_mat, paste0(output_dir,"/normal.txt"), sep="\t", quote=FALSE)
write.table(data.frame(cell_names), paste0(output_dir,"/barcodes.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("  Wrote tumor.txt, normal.txt, barcodes.tsv\n")

# ── SNP sparse matrices (vectorized) ─────────────────────────────────────────
set.seed(42)
snps_per_bin <- 2
n_snps <- n_bins * snps_per_bin
cat("  Building", n_snps, "SNPs x", n_cells, "cells sparse matrices...\n")

snp_bin_idx <- rep(seq_len(n_bins), each=snps_per_bin)  # which bin each SNP belongs to

alt_mat_full <- as.matrix(alt_bin)[, , drop=FALSE]  # n_cells x n_bins
ref_mat_full <- as.matrix(ref_bin)[, , drop=FALSE]

# Expand bin counts to SNP level by repeating columns
alt_snp <- t(alt_mat_full[, snp_bin_idx, drop=FALSE])  # n_snps x n_cells
ref_snp <- t(ref_mat_full[, snp_bin_idx, drop=FALSE])

# Distribute counts evenly across SNPs per bin
alt_snp <- round(alt_snp / snps_per_bin)
ref_snp <- round(ref_snp / snps_per_bin)

writeMM(Matrix(alt_snp, sparse=TRUE), paste0(output_dir,"/alt_all.mtx"))
writeMM(Matrix(ref_snp, sparse=TRUE), paste0(output_dir,"/ref_all.mtx"))
cat("  Wrote alt_all.mtx and ref_all.mtx\n")

# ── VCF ───────────────────────────────────────────────────────────────────────
snp_positions <- as.integer(unlist(lapply(seq_len(n_bins), function(bi) {
    seq(bins$start[bi]+1, bins$end[bi]-1, length.out=snps_per_bin)
})))
vcf_df <- data.frame(
    V1=rep(bins$chrom, each=snps_per_bin), V2=snp_positions,
    V3=".", V4="A", V5="T", V6=".", V7="PASS", V8=".",
    stringsAsFactors=FALSE
)
write.table(vcf_df, paste0(output_dir,"/var_all.vcf"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("  Wrote var_all.vcf\n")

if (length(normal_cells) > 0)
    writeLines(normal_cells, paste0(output_dir,"/normal_cells.txt"))

cat("Done. Output in", output_dir, "\n")
