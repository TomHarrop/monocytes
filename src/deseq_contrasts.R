#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(DESeq2)
library(data.table)

sample_data_file <- snakemake@input[["factor_data"]]
dds_file <-  snakemake@input[["dds"]]
ncpus <- snakemake@threads
lfc_thresh <- snakemake@params[["lfc_threshold"]]
alpha <- snakemake@params[["alpha"]]

# dev
# sample_data_file <- "data/factor_data.csv"
# dds_file <- "output/030_deseq/dds.Rds"
# ncpus <- 8
# lfc_thresh <- log(1.5, 2)
# alpha <- 0.1

# register parallel processors
BiocParallel::register(BiocParallel::MulticoreParam(ncpus))

# sample_data
sample_data <- fread(sample_data_file)

# load deseq object
dds_all <- readRDS(dds_file)

dds_with_sd <- DESeqDataSetFromMatrix(
    countData = counts(dds_all),
    colData = data.frame(sample_data, row.names = "Sample"),
    design = ~ Donor + Treatment
)
    
    
# at least 3 samples with a count of 10 or higher
# chucks out a handful more genes
keep_genes <- rowSums(counts(dds_with_sd) >= 10) >= 3
dds_filtered <- dds_with_sd[keep_genes,]

# run deseq
dds <- DESeq(dds_filtered, parallel = TRUE)

# make contrasts
contrasts <- list(
    c("Treatment", "Untreated", "Intact"),
    c("Treatment", "Untreated", "Lysed"),
    c("Treatment", "Lysed", "Intact"))

contrasts <- c(contrasts,
               lapply(contrasts, function(x) x[c(1, 3, 2)]))


# extract results
all_res <- lapply(
    contrasts, function(x)
    results(object = dds,
            contrast = x,
            lfcThreshold = lfc_thresh,
            alpha = alpha, tidy = TRUE))

names(all_res) <- sapply(contrasts, paste, collapse = ".")
res_dt <- rbindlist(all_res, idcol = "contrast")
setnames(res_dt, "row", "gene")
setorder(res_dt, gene, contrast)

# get counts
norm_counts <- data.table(as.data.frame(counts(dds, normalized = TRUE)),
                       keep.rownames = "gene")

# subset counts by comparison
contrast_gene_lists <- lapply(all_res, function(x)
    norm_counts[gene %in% data.table(x)[padj < alpha, unique(row)]])
contrast_genes <- rbindlist(contrast_gene_lists, idcol = "contrast")
setorder(contrast_genes, contrast, gene)

# write output
fwrite(res_dt, snakemake@output[["wald_results"]])
fwrite(norm_counts, snakemake@output[["counts"]])
fwrite(contrast_genes, snakemake@output[["contrast_genes"]])

# log
sessionInfo()
