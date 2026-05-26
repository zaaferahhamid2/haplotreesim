# Evaluate Alleloscope output against HaploTreeSim ground truth
# Usage: Rscript scripts/evaluate_alleloscope.R \
#          --dataset-dir examples/whole_genome_simulation \
#          --output-dir examples/whole_genome_simulation/alleloscope_output \
#          --metrics-out examples/whole_genome_simulation/alleloscope_metrics.json

args <- commandArgs(trailingOnly=TRUE)
dataset_dir  <- args[which(args=="--dataset-dir")  + 1]
output_dir   <- args[which(args=="--output-dir")   + 1]
metrics_out  <- args[which(args=="--metrics-out")  + 1]
if (is.na(metrics_out)) metrics_out <- paste0(output_dir, "/metrics.json")

cat("Loading ground truth from", dataset_dir, "\n")
truth   <- read.table(paste0(dataset_dir,"/truth_cell_hscn_segments.tsv"), sep="\t", header=TRUE)
cells   <- read.table(paste0(dataset_dir,"/cells.tsv"), sep="\t", header=TRUE)
segments <- read.table(paste0(dataset_dir,"/segments.tsv"), sep="\t", header=TRUE)

n_cells   <- nrow(cells)
n_segs    <- nrow(segments)
cell_names <- cells$cell
clone_labels <- cells$clone_assignment

cat("  Cells:", n_cells, "| Segments:", n_segs, "\n")

# True CN matrices
cn_A_true <- matrix(0L, nrow=n_cells, ncol=n_segs)
cn_B_true <- matrix(0L, nrow=n_cells, ncol=n_segs)
for (i in seq_len(nrow(truth))) {
    ci <- truth$cell_index[i] + 1
    si <- truth$segment[i] + 1
    cn_A_true[ci, si] <- truth$cn_A[i]
    cn_B_true[ci, si] <- truth$cn_B[i]
}

# Load Alleloscope genotype_values (rho=total CN ratio, theta=BAF)
cat("Loading Alleloscope output...\n")
gv_file <- paste0(output_dir, "/genotype_values.tsv")
if (!file.exists(gv_file)) stop("genotype_values.tsv not found")
gv <- read.table(gv_file, sep="\t", header=TRUE)
cat("  genotype_values:", nrow(gv), "cells x", ncol(gv), "columns\n")

# Extract rho and theta per cell per region
regions <- unique(gsub("rho_|theta_|h1_|h2_", "", colnames(gv)))
# R converts ":" to "." in column names
regions <- gsub("\\.", ":", regions)
regions <- regions[grepl(":", regions)]
cat("  Regions:", paste(regions, collapse=", "), "\n")

# Since we have 1 region (whole chr1), map to all segments
# rho ~ total CN / ploidy, theta ~ BAF
rho_col   <- paste0("rho_", gsub(":", ".", regions[1]))
theta_col <- paste0("theta_", gsub(":", ".", regions[1]))

if (!rho_col %in% colnames(gv)) {
    cat("WARNING: rho column not found, using 1.0\n")
    rho_vals <- rep(1.0, n_cells)
} else {
    rho_vals <- gv[match(cell_names, rownames(gv)), rho_col]
    rho_vals[is.na(rho_vals)] <- 1.0
}

if (!theta_col %in% colnames(gv)) {
    cat("WARNING: theta column not found, using 0.5\n")
    theta_vals <- rep(0.5, n_cells)
} else {
    theta_vals <- gv[match(cell_names, rownames(gv)), theta_col]
    theta_vals[is.na(theta_vals)] <- 0.5
}

# Convert rho/theta to haplotype CN
# rho = total_CN / 2 (diploid baseline), theta = cn_A / total_CN
# So: total_CN = rho * 2, cn_A = theta * total_CN, cn_B = (1-theta) * total_CN
total_cn_pred <- pmax(round(rho_vals * 2), 0)
cn_A_pred <- matrix(0L, nrow=n_cells, ncol=n_segs)
cn_B_pred <- matrix(0L, nrow=n_cells, ncol=n_segs)
for (si in seq_len(n_segs)) {
    cn_A_pred[, si] <- pmax(round(theta_vals * total_cn_pred), 0L)
    cn_B_pred[, si] <- pmax(total_cn_pred - cn_A_pred[, si], 0L)
}

# ── Metrics ───────────────────────────────────────────────────────────────────
# HSCN swap-invariant error
hscn_err <- mean(pmin(
    abs(cn_A_pred - cn_A_true) + abs(cn_B_pred - cn_B_true),
    abs(cn_A_pred - cn_B_true) + abs(cn_B_pred - cn_A_true)
))

# LOH F1
true_loh <- (pmin(cn_A_true, cn_B_true) == 0) & ((cn_A_true + cn_B_true) >= 1)
pred_loh <- (pmin(cn_A_pred, cn_B_pred) == 0) & ((cn_A_pred + cn_B_pred) >= 1)
tp <- sum(true_loh & pred_loh)
fp <- sum(!true_loh & pred_loh)
fn <- sum(true_loh & !pred_loh)
loh_precision <- if (tp+fp > 0) tp/(tp+fp) else NA
loh_recall    <- if (tp+fn > 0) tp/(tp+fn) else NA
loh_f1        <- if (!is.na(loh_precision) && !is.na(loh_recall) && (loh_precision+loh_recall)>0)
    2*loh_precision*loh_recall/(loh_precision+loh_recall) else NA

# Clone ARI (simple: use total CN profile for clustering)
tcn_pred <- cn_A_pred + cn_B_pred
if (requireNamespace("mclust", quietly=TRUE)) {
    ari <- mclust::adjustedRandIndex(clone_labels, 
           kmeans(tcn_pred, centers=length(unique(clone_labels)), nstart=10)$cluster)
} else {
    ari <- NA
    cat("  Note: install mclust for ARI\n")
}

cat("\n==================================================\n")
cat("ALLELOSCOPE EVALUATION RESULTS\n")
cat("==================================================\n")
cat(sprintf("  HSCN Error:  %.4f  (0=perfect)\n", hscn_err))
cat(sprintf("  LOH F1:      %.4f\n", ifelse(is.na(loh_f1), -1, loh_f1)))
cat(sprintf("  Clone ARI:   %.4f  (1=perfect)\n", ifelse(is.na(ari), -1, ari)))
cat("==================================================\n")
cat("\nTrue CN sample (cell 1):", paste(cn_A_true[1,], cn_B_true[1,], sep=","), "\n")
cat("Pred CN sample (cell 1):", paste(cn_A_pred[1,], cn_B_pred[1,], sep=","), "\n")

# Save metrics
metrics <- list(
    hscn_error=hscn_err,
    loh_precision=loh_precision,
    loh_recall=loh_recall,
    loh_f1=loh_f1,
    clone_ari=ari,
    n_cells=n_cells,
    n_segments=n_segs,
    n_regions_alleloscope=length(regions)
)
json_str <- paste0('{\n', paste(sapply(names(metrics), function(k) {
    v <- metrics[[k]]
    if (is.na(v)) sprintf('  "%s": null', k)
    else if (is.numeric(v)) sprintf('  "%s": %.6f', k, v)
    else sprintf('  "%s": %s', k, v)
}), collapse=',\n'), '\n}')
writeLines(json_str, metrics_out)
cat("Wrote metrics to", metrics_out, "\n")
