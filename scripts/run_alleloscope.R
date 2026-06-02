# Run Alleloscope on converted HaploTreeSim input
# Usage: Rscript scripts/run_alleloscope.R \
#          --input-dir examples/whole_genome_simulation/alleloscope_input \
#          --output-dir examples/whole_genome_simulation/alleloscope_output

library(Alleloscope)
library(Matrix)

args <- commandArgs(trailingOnly=TRUE)
input_dir  <- args[which(args=="--input-dir")  + 1]
output_dir <- args[which(args=="--output-dir") + 1]

dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(paste0(output_dir,"/rds"), showWarnings=FALSE)
dir.create(paste0(output_dir,"/rds/EMresults"), showWarnings=FALSE)
dir.create(paste0(output_dir,"/plots/EMresults"), showWarnings=FALSE)
dir.create(paste0(output_dir,"/plots"), showWarnings=FALSE)

cat("Loading Alleloscope inputs from", input_dir, "\n")
barcodes   <- read.table(paste0(input_dir,"/barcodes.tsv"), stringsAsFactors=FALSE, header=FALSE)
alt_all    <- readMM(paste0(input_dir,"/alt_all.mtx"))
ref_all    <- readMM(paste0(input_dir,"/ref_all.mtx"))
var_all    <- read.table(paste0(input_dir,"/var_all.vcf"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
raw_counts <- read.table(paste0(input_dir,"/tumor.txt"),  sep="\t", header=TRUE, row.names=1)
ref_counts <- read.table(paste0(input_dir,"/normal.txt"), sep="\t", header=TRUE, row.names=1)

# Use normal cells if available
normal_cells_file <- paste0(input_dir,"/normal_cells.txt")
normal_cells <- if (file.exists(normal_cells_file)) readLines(normal_cells_file) else NULL
cat("  Normal cells available:", length(normal_cells), "\n")

# Chromosome sizes from data
chroms <- unique(var_all$V1)
chrom_nums <- gsub("chr","", chroms)
size <- data.frame(
    V1=chrom_nums,
    V2=sapply(chrom_nums, function(c) {
        rows <- raw_counts[grepl(paste0("chr",c,"-"), rownames(raw_counts)),]
        if(nrow(rows)==0) return(248956422L)
        parts <- strsplit(rownames(rows)[nrow(rows)], "-")[[1]]
        as.integer(parts[3])
    }),
    stringsAsFactors=FALSE
)

data(centromere.GRCh38)
data(telomere.GRCh38)

cat("Running Alleloscope...\n")
obj <- Alleloscope:::Createobj(
    alt_all=alt_all, ref_all=ref_all, var_all=var_all,
    samplename="HaploTreeSim", genome_assembly="GRCh38",
    dir_path=paste0(output_dir,"/"),
    barcodes=barcodes, size=size, assay="scDNAseq"
)

obj <- Alleloscope:::Matrix_filter(obj, cell_filter=1, SNP_filter=2,
    min_vaf=0.05, max_vaf=0.95,
    centro=centromere.GRCh38, telo=telomere.GRCh38)
cat("  Filtered. Cells:", ncol(obj$total_all), "SNPs:", nrow(obj$total_all), "\n")

obj <- Alleloscope:::Segmentation(Obj_filtered=obj,
    raw_counts=raw_counts, ref_counts=ref_counts,
    rds_path=paste0(output_dir,"/rds/genotype_"))
cat("  Segments found:", nrow(obj$seg_table), "\n")

obj <- Alleloscope:::Segments_filter(Obj_filtered=obj, nSNP=2, len=100000)
cat("  Segments after filter:", nrow(obj$seg_table_filtered), "\n")

if (is.null(obj$seg_table_filtered) || nrow(obj$seg_table_filtered) == 0) {
    cat("  WARNING: No segments found. Forcing whole chromosome.\n")
    obj$seg_table_filtered <- obj$seg_table
    if (!is.null(obj$seg_table_filtered) && nrow(obj$seg_table_filtered) > 0)
        obj$seg_table_filtered$chrr <- paste0(obj$seg_table_filtered$chr, ":0")
}

obj <- Alleloscope:::Est_regions(Obj_filtered=obj,
    rds_path=paste0(output_dir,"/rds/genotype_"),
    max_nSNP=30000, min_cell=1, min_snp=0)
cat("  EM estimation done\n")

obj <- Alleloscope:::Select_normal(Obj_filtered=obj,
    raw_counts=raw_counts, cell_nclust=5, plot_theta=FALSE,
    pre_sel=!is.null(normal_cells),
    cell_type=if(!is.null(normal_cells)) {
        ct <- rep("unknown", ncol(obj$total_all))
        names(ct) <- colnames(obj$total_all)
        ct[names(ct) %in% normal_cells] <- "normal"
        ct
    } else NULL)
cat("  Select_normal done\n")

# Set ref to the normal region for Genotype_value
if (is.null(obj$ref)) obj$ref <- obj$seg_table_filtered$chrr[1]
obj <- Alleloscope:::Genotype_value(Obj_filtered=obj, type="tumor", refr=FALSE,
    raw_counts=raw_counts, ref_counts=ref_counts)
cat("  Genotype_value done\n")

obj$dir_path <- paste0(output_dir, "/")
obj <- Alleloscope:::Genotype(Obj_filtered=obj,
    rds_path=NULL, plot_path=NULL)
cat("  Genotype done\n")

saveRDS(obj, paste0(output_dir,"/rds/obj_final.rds"))

# Save genotypes as TSV for evaluator
geno <- obj$genotypes
if (!is.null(geno)) {
    write.table(geno, paste0(output_dir,"/genotypes.tsv"), sep="\t", quote=FALSE)
    cat("  Wrote genotypes.tsv:", nrow(geno), "cells x", ncol(geno), "regions\n")
}

gv <- obj$genotype_values
if (!is.null(gv)) {
    write.table(gv, paste0(output_dir,"/genotype_values.tsv"), sep="\t", quote=FALSE)
    cat("  Wrote genotype_values.tsv\n")
}

cat("Done. Output in", output_dir, "\n")
