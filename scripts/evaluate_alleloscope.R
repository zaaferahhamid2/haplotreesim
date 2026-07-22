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
regions <- sub("\\.", ":", regions)  # only first "." is the chrom:pos separator; later periods may be decimals in sci notation
regions <- regions[grepl(":", regions)]
cat("  Regions:", paste(regions, collapse=", "), "\n")

# Alleloscope regions are chromosome-level in the current runner (for example
# "1:0"). Map each truth segment to the matching chromosome region instead of
# applying the first region to every segment.
region_chrom <- paste0("chr", sub(":.*", "", regions))
region_by_chrom <- setNames(regions, region_chrom)
cell_match <- match(cell_names, rownames(gv))

cn_A_pred <- matrix(0L, nrow=n_cells, ncol=n_segs)
cn_B_pred <- matrix(0L, nrow=n_cells, ncol=n_segs)
for (si in seq_len(n_segs)) {
    seg_chrom <- as.character(segments$chrom[si])
    region <- region_by_chrom[[seg_chrom]]
    if (is.null(region) || is.na(region)) {
        rho_vals <- rep(1.0, n_cells)
        theta_vals <- rep(0.5, n_cells)
    } else {
        rho_col <- paste0("rho_", gsub(":", ".", region))
        theta_col <- paste0("theta_", gsub(":", ".", region))

        if (!rho_col %in% colnames(gv)) {
            rho_vals <- rep(1.0, n_cells)
        } else {
            rho_vals <- gv[cell_match, rho_col]
            rho_vals[is.na(rho_vals)] <- 1.0
        }

        if (!theta_col %in% colnames(gv)) {
            theta_vals <- rep(0.5, n_cells)
        } else {
            theta_vals <- gv[cell_match, theta_col]
            theta_vals[is.na(theta_vals)] <- 0.5
        }
    }

    # rho = total_CN / 2 (diploid baseline), theta = cn_A / total_CN.
    total_cn_pred <- pmax(round(rho_vals * 2), 0)
    cn_A_pred[, si] <- pmax(round(theta_vals * total_cn_pred), 0L)
    cn_B_pred[, si] <- pmax(total_cn_pred - cn_A_pred[, si], 0L)
}

pred_cell_names <- rownames(gv)
common_cells <- intersect(cell_names, pred_cell_names)
cat('  Cells with predictions:', length(common_cells), '\n')
cell_mask_retained <- cell_names %in% common_cells
n_cells_retained <- sum(cell_mask_retained)
retention_rate <- n_cells_retained / n_cells

retention_by_clone <- aggregate(
    cell_mask_retained,
    by=list(truth_clone=cells$clone),
    FUN=function(x) c(n_total=length(x), n_retained=sum(x))
)
retention_by_clone <- data.frame(
    truth_clone=retention_by_clone$truth_clone,
    n_total=retention_by_clone$x[, "n_total"],
    n_retained=retention_by_clone$x[, "n_retained"],
    stringsAsFactors=FALSE
)
retention_by_clone$retention_rate <- ifelse(
    retention_by_clone$n_total > 0,
    retention_by_clone$n_retained / retention_by_clone$n_total,
    NA
)

# ── Metrics ───────────────────────────────────────────────────────────────────
compute_scope_metrics <- function(mask) {
    n_scope <- sum(mask)
    if (n_scope == 0) {
        return(list(
            n_cells=0L,
            hscn_error=NA,
            loh_precision=NA,
            loh_recall=NA,
            loh_f1=NA,
            loh_n_true=0L,
            loh_n_pred=0L,
            loh_tp=0L,
            loh_fp=0L,
            loh_fn=0L,
            clone_ari=NA
        ))
    }

    pred_A <- cn_A_pred[mask, , drop=FALSE]
    pred_B <- cn_B_pred[mask, , drop=FALSE]
    true_A <- cn_A_true[mask, , drop=FALSE]
    true_B <- cn_B_true[mask, , drop=FALSE]
    labels <- clone_labels[mask]

    hscn_err <- mean(pmin(
        abs(pred_A - true_A) + abs(pred_B - true_B),
        abs(pred_A - true_B) + abs(pred_B - true_A)
    ))

    true_loh <- (pmin(true_A, true_B) == 0) & ((true_A + true_B) >= 1)
    pred_loh <- (pmin(pred_A, pred_B) == 0) & ((pred_A + pred_B) >= 1)
    tp <- sum(true_loh & pred_loh)
    fp <- sum(!true_loh & pred_loh)
    fn <- sum(true_loh & !pred_loh)
    n_true_loh <- sum(true_loh)
    n_pred_loh <- sum(pred_loh)
    if (n_true_loh == 0) {
        loh_precision <- NA
        loh_recall <- NA
        loh_f1 <- NA
    } else {
        loh_precision <- if (tp + fp > 0) tp / (tp + fp) else 0.0
        loh_recall <- tp / (tp + fn)
        loh_f1 <- if (loh_precision + loh_recall > 0)
            2 * loh_precision * loh_recall / (loh_precision + loh_recall) else 0.0
    }

    tcn_pred <- pred_A + pred_B
    if (requireNamespace('mclust', quietly=TRUE)) {
        if (nrow(unique(tcn_pred)) < 2) {
            pred_cluster <- rep(1L, nrow(tcn_pred))
        } else {
            n_centers <- min(length(unique(labels)), nrow(unique(tcn_pred)))
            set.seed(42)
            pred_cluster <- kmeans(tcn_pred, centers=n_centers, nstart=10)$cluster
        }
        ari <- mclust::adjustedRandIndex(labels, pred_cluster)
    } else {
        ari <- NA
        cat('  Note: mclust unavailable for ARI\n')
    }

    list(
        n_cells=as.integer(n_scope),
        hscn_error=hscn_err,
        loh_precision=loh_precision,
        loh_recall=loh_recall,
        loh_f1=loh_f1,
        loh_n_true=as.integer(n_true_loh),
        loh_n_pred=as.integer(n_pred_loh),
        loh_tp=as.integer(tp),
        loh_fp=as.integer(fp),
        loh_fn=as.integer(fn),
        clone_ari=ari
    )
}

all_cells_metrics <- compute_scope_metrics(rep(TRUE, n_cells))
retained_cells_metrics <- compute_scope_metrics(cell_mask_retained)

cat("\n==================================================\n")
cat("ALLELOSCOPE EVALUATION RESULTS\n")
cat("==================================================\n")
cat(sprintf("  Retention:   %d/%d cells (%.4f)\n",
            n_cells_retained, n_cells, retention_rate))
cat(sprintf("  All cells HSCN Error:       %.4f  (0=perfect)\n",
            all_cells_metrics$hscn_error))
cat(sprintf("  All cells LOH F1:           %s  (true=%d pred=%d)\n",
            ifelse(is.na(all_cells_metrics$loh_f1), "NA", sprintf("%.4f", all_cells_metrics$loh_f1)),
            all_cells_metrics$loh_n_true, all_cells_metrics$loh_n_pred))
cat(sprintf("  All cells Clone ARI:        %s  (1=perfect)\n",
            ifelse(is.na(all_cells_metrics$clone_ari), "NA", sprintf("%.4f", all_cells_metrics$clone_ari))))
cat(sprintf("  Retained cells HSCN Error:  %.4f  (0=perfect)\n",
            retained_cells_metrics$hscn_error))
cat(sprintf("  Retained cells LOH F1:      %s  (true=%d pred=%d)\n",
            ifelse(is.na(retained_cells_metrics$loh_f1), "NA", sprintf("%.4f", retained_cells_metrics$loh_f1)),
            retained_cells_metrics$loh_n_true, retained_cells_metrics$loh_n_pred))
cat(sprintf("  Retained cells Clone ARI:   %s  (1=perfect)\n",
            ifelse(is.na(retained_cells_metrics$clone_ari), "NA", sprintf("%.4f", retained_cells_metrics$clone_ari))))
cat("==================================================\n")
cat("\nTrue CN sample (cell 1):", paste(cn_A_true[1,], cn_B_true[1,], sep=","), "\n")
cat("Pred CN sample (cell 1):", paste(cn_A_pred[1,], cn_B_pred[1,], sep=","), "\n")

# Save metrics
metrics <- list(
    all_cells=all_cells_metrics,
    retained_cells=retained_cells_metrics,
    retention_by_truth_clone=retention_by_clone,
    retention_rate=retention_rate,
    n_cells=n_cells,
    n_cells_retained=n_cells_retained,
    n_segments=n_segs,
    n_regions_alleloscope=length(regions),
    # Legacy top-level keys: HSCN/LOH are all-cell metrics; clone_ari keeps
    # the previous retained-cell behavior.
    hscn_error=all_cells_metrics$hscn_error,
    loh_precision=all_cells_metrics$loh_precision,
    loh_recall=all_cells_metrics$loh_recall,
    loh_f1=all_cells_metrics$loh_f1,
    loh_n_true=all_cells_metrics$loh_n_true,
    loh_n_pred=all_cells_metrics$loh_n_pred,
    loh_tp=all_cells_metrics$loh_tp,
    loh_fp=all_cells_metrics$loh_fp,
    loh_fn=all_cells_metrics$loh_fn,
    clone_ari=retained_cells_metrics$clone_ari
)
json_str <- jsonlite::toJSON(
    metrics,
    pretty=TRUE,
    auto_unbox=TRUE,
    dataframe="rows",
    na="null"
)
writeLines(json_str, metrics_out)
cat("Wrote metrics to", metrics_out, "\n")
