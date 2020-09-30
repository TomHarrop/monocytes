#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)
library(DESeq2)
library(rtracklayer)
library(tximport)


file_list <- snakemake@input[["quant_files"]]
gff <- snakemake@input[["gff"]]
dds_file <- snakemake@output[[1]]

# dev
 # gff <- "data/ref/GCF_000001405.39_GRCh38.p13_genomic.gff"
 # file_list <- list.files("output/020_salmon",
 #                         full.names = TRUE,
 #                         recursive = TRUE,
 #                         pattern = "quant.sf")


# read gff using GenomicRanges
gr <- import.gff3(gff,
                  feature.type = c("exon",
                                   "CDS",
                                   "mRNA",
                                   "gene",
                                   "transcript"))

# extract a data.frame of tx to gene
mrnas <- gr[gr$type %in% c("mRNA", "transcript"),]
mrna_dt <- as.data.table(mcols(mrnas))
tx2gene <- data.frame(mrna_dt[, .(TXNAME = Name,
                                  GENEID = as.character(gene))])

# get file list
names(file_list) <- gsub(".*/", "", dirname(file_list))

# import quant files
txi <- tximport(file_list, type = "salmon", tx2gene = tx2gene)

# which stuff isn't in tx2gene?
read_list <- lapply(file_list, fread)
reads <- rbindlist(read_list, idcol = "sample")
missing_genes <- reads[!Name %in% tx2gene$TXNAME, unique(Name)]

# generate col_data
col_data <- data.table(samplename = names(file_list))
col_data[, group := samplename]

# generate DESeq object
dds <- DESeqDataSetFromTximport(
  txi,
  colData = data.frame(col_data, row.names = 'samplename'),
  design = ~ group
)

saveRDS(dds, dds_file)

# log
sessionInfo()
