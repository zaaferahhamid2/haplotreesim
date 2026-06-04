library(Matrix)

args <- commandArgs(trailingOnly=TRUE)

arg_value <- function(name) {
    idx <- which(args == name)
    if (length(idx) == 0 || idx[1] == length(args)) {
        stop("Missing required argument ", name, call.=FALSE)
    }
    args[idx[1] + 1]
}

input_dir  <- arg_value("--input-dir")
output_dir <- arg_value("--output-dir")

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings=FALSE)
dir.create(file.path(output_dir, "rds"), showWarnings=FALSE)

required_files <- c(
    "cells.tsv",
    "bins.tsv",
    "readcounts.tsv",
    "snps.tsv",
    "snp_allele_alt.mtx",
    "snp_allele_ref.mtx"
)
missing_files <- required_files[!file.exists(file.path(input_dir, required_files))]
if (length(missing_files) > 0) {
    stop(
        "Alleloscope conversion requires SNP-level simulator outputs. Missing: ",
        paste(missing_files, collapse=", "),
        call.=FALSE
    )
}

cat("Loading dataset...\n")
cells <- read.table(file.path(input_dir, "cells.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
bins  <- read.table(file.path(input_dir, "bins.tsv"),  sep="\t", header=TRUE, stringsAsFactors=FALSE)
rc    <- read.table(file.path(input_dir, "readcounts.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE)
snps  <- read.table(file.path(input_dir, "snps.tsv"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

cat("  Reading SNP-level sparse allele matrices...\n")
alt_snp <- readMM(file.path(input_dir, "snp_allele_alt.mtx"))
ref_snp <- readMM(file.path(input_dir, "snp_allele_ref.mtx"))

n_cells    <- nrow(cells)
n_bins     <- nrow(bins)
cell_names <- cells$cell
normal_cells <- cells$cell[cells$clone_assignment == 0]
n_snps     <- nrow(snps)

if (!all(c("chrom", "position", "bin") %in% colnames(snps))) {
    stop("snps.tsv must contain chrom, position, and bin columns", call.=FALSE)
}
if (!all(dim(alt_snp) == dim(ref_snp))) {
    stop("snp_allele_alt.mtx and snp_allele_ref.mtx have different dimensions", call.=FALSE)
}
if (nrow(alt_snp) != n_snps) {
    stop(
        "SNP matrix row count (", nrow(alt_snp), ") does not match snps.tsv rows (", n_snps, ")",
        call.=FALSE
    )
}
if (ncol(alt_snp) != n_cells) {
    stop(
        "SNP matrix column count (", ncol(alt_snp), ") does not match cells.tsv rows (", n_cells, ")",
        call.=FALSE
    )
}
if (any(is.na(snps$position))) {
    stop("snps.tsv contains missing SNP positions", call.=FALSE)
}
if (any(alt_snp@x < 0) || any(ref_snp@x < 0)) {
    stop("SNP-level allele matrices contain negative counts", call.=FALSE)
}

cat("  Cells:", n_cells, "| Bins:", n_bins, "| SNPs:", n_snps, "| Normal:", length(normal_cells), "\n")

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

# ── SNP sparse matrices ──────────────────────────────────────────────────────
writeMM(alt_snp, paste0(output_dir,"/alt_all.mtx"))
writeMM(ref_snp, paste0(output_dir,"/ref_all.mtx"))
cat("  Wrote real SNP-level alt_all.mtx and ref_all.mtx\n")

# ── VCF ───────────────────────────────────────────────────────────────────────
vcf_df <- data.frame(
    V1=as.character(snps$chrom), V2=as.integer(snps$position),
    V3=".", V4="A", V5="T", V6=".", V7="PASS", V8=".",
    stringsAsFactors=FALSE
)
write.table(vcf_df, paste0(output_dir,"/var_all.vcf"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("  Wrote var_all.vcf from snps.tsv\n")

if (length(normal_cells) > 0)
    writeLines(normal_cells, paste0(output_dir,"/normal_cells.txt"))

cat("Done. Output in", output_dir, "\n")
