# Alleloscope integration script (Week 16)
# Status: Partial - segmentation step finds only 1 segment due to flat RDR
# (no matched normal sample available from simulator)
# Alleloscope runs successfully on its own sample data (SNU601)

library(Alleloscope)
library(Matrix)

setwd("~/Documents/haplotreesim")
dir_path <- "alleloscope_output/"
dir.create(dir_path, showWarnings=FALSE)
dir.create(paste0(dir_path, "rds"), showWarnings=FALSE)
dir.create(paste0(dir_path, "plots"), showWarnings=FALSE)

cat("Loading simulated data...\n")
rc  <- read.table("seacon_output/readcounts.tsv", sep="\t", header=TRUE, row.names=1)
baf <- read.table("seacon_output/precomputed_baf.tsv", sep="\t", header=TRUE)

n_cells  <- nrow(rc)
n_bins   <- ncol(rc)
bin_size <- 500000
cell_names <- rownames(rc)
cat("  Cells:", n_cells, "| Bins:", n_bins, "\n")

bin_starts <- formatC((0:(n_bins-1)) * bin_size + 1, format="d")
bin_ends   <- formatC((1:n_bins) * bin_size, format="d")
row_names  <- paste0("chr1-", bin_starts, "-", bin_ends)

tumor_mat <- as.data.frame(t(rc))
rownames(tumor_mat) <- row_names
colnames(tumor_mat) <- cell_names
tumor_mat <- tumor_mat[, , drop=FALSE]

normal_mean <- round(mean(colMeans(tumor_mat)))
set.seed(42)
normal_vec <- round(rnorm(n_bins, normal_mean, normal_mean * 0.05))
normal_mat <- data.frame(normal=pmax(normal_vec, 1), row.names=row_names)

write.table(tumor_mat,  paste0(dir_path, "tumor.txt"),  sep="\t", quote=FALSE)
write.table(normal_mat, paste0(dir_path, "normal_diploid.txt"), sep="\t", quote=FALSE)
cat("  Wrote tumor.txt and normal_diploid.txt\n")

set.seed(42)
snps_per_seg <- 200
segments <- unique(baf[, c("chrom","start","end")])
n_segs <- nrow(segments)
n_snps <- n_segs * snps_per_seg
cat("  Segments:", n_segs, "| Total SNPs:", n_snps, "\n")

vcf_rows <- data.frame(
  V1=rep("chr1", n_snps),
  V2=as.integer(unlist(lapply(1:n_segs, function(si) {
    seg <- segments[si,]
    as.integer(seq(seg$start+1, seg$end-1, length.out=snps_per_seg))
  }))),
  V3=".", V4="A", V5="T", V6=".", V7="PASS", V8=".",
  stringsAsFactors=FALSE
)

snp_seg_idx <- rep(1:n_segs, each=snps_per_seg)
baf_mat <- matrix(0.5, nrow=n_snps, ncol=n_cells)
for (si in 1:n_segs) {
  seg <- segments[si,]
  seg_baf <- baf[baf$start == seg$start,]
  for (ci in 1:n_cells) {
    cell_row <- seg_baf[seg_baf$cell == cell_names[ci],]
    if (nrow(cell_row) > 0)
      baf_mat[snp_seg_idx == si, ci] <- cell_row$BAF[1]
  }
}

total_mat <- matrix(rpois(n_snps * n_cells, 10), nrow=n_snps, ncol=n_cells)
total_mat <- pmax(total_mat, 1L)
alt_mat_full <- matrix(rbinom(n_snps * n_cells, as.vector(total_mat),
                               as.vector(baf_mat)), nrow=n_snps, ncol=n_cells)
ref_mat_full <- total_mat - alt_mat_full

alt_sparse <- Matrix(alt_mat_full, sparse=TRUE, doDiag=FALSE)
ref_sparse <- Matrix(ref_mat_full, sparse=TRUE, doDiag=FALSE)
writeMM(alt_sparse, paste0(dir_path, "alt_all.mtx"))
writeMM(ref_sparse, paste0(dir_path, "ref_all.mtx"))

write.table(vcf_rows, paste0(dir_path, "var_all.vcf"), sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(data.frame(cell_names), paste0(dir_path, "barcodes.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("  Wrote all input files\n")

cat("\nRunning Alleloscope step by step...\n")
data(centromere.GRCh38)
data(telomere.GRCh38)

size     <- data.frame(V1="1", V2=248956422, stringsAsFactors=FALSE)
barcodes <- read.table(paste0(dir_path,"barcodes.tsv"), stringsAsFactors=FALSE, header=FALSE)
alt_all  <- readMM(paste0(dir_path,"alt_all.mtx"))
ref_all  <- readMM(paste0(dir_path,"ref_all.mtx"))
var_all  <- read.table(paste0(dir_path,"var_all.vcf"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
raw_counts <- read.table(paste0(dir_path,"tumor.txt"), sep="\t", header=TRUE, row.names=1)
ref_counts <- read.table(paste0(dir_path,"normal_diploid.txt"), sep="\t", header=TRUE, row.names=1)

obj <- Alleloscope:::Createobj(alt_all=alt_all, ref_all=ref_all, var_all=var_all,
    samplename="HaploTreeSim", genome_assembly="GRCh38", dir_path=dir_path,
    barcodes=barcodes, size=size, assay="scDNAseq")

obj <- Alleloscope:::Matrix_filter(obj, cell_filter=1, SNP_filter=5,
    min_vaf=0.05, max_vaf=0.95, centro=centromere.GRCh38, telo=telomere.GRCh38)
cat("Filtered. Cells:", ncol(obj$total_all), "SNPs:", nrow(obj$total_all), "\n")

obj <- Alleloscope:::Segmentation(Obj_filtered=obj, raw_counts=raw_counts,
    ref_counts=ref_counts, rds_path=paste0(dir_path,"rds/"))
cat("Segments found:", nrow(obj$seg_table), "\n")

obj <- Alleloscope:::Segments_filter(Obj_filtered=obj, nSNP=5, len=100000)
cat("Segments after filter:", nrow(obj$seg_table_filtered), "\n")

# Force continue with 1 segment if needed
if (is.null(obj$seg_table_filtered) || nrow(obj$seg_table_filtered) == 0) {
  cat("WARNING: No segments after filter. Forcing whole chromosome as 1 segment.\n")
  obj$seg_table_filtered <- obj$seg_table
  obj$seg_table_filtered$chrr <- "1:0"
}

obj <- Alleloscope:::Est_regions(Obj_filtered=obj,
    rds_path=paste0(dir_path,"rds/"), max_nSNP=30000, min_cell=1, min_snp=0)
cat("EM estimation done\n")

obj <- Alleloscope:::Select_normal(Obj_filtered=obj, raw_counts=raw_counts,
    cell_nclust=5, plot_theta=FALSE)
cat("Select_normal done\n")

saveRDS(obj, paste0(dir_path,"rds/obj_final.rds"))
cat("\nDone. Output in", dir_path, "\n")
