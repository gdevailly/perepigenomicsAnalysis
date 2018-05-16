library(dplyr)
library(readr)
library(ggplot2)
setwd("/groups2/joshi_grp/guillaume/cascade/")

tss <- read_tsv("../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.bed", col_names = FALSE)
tes <- read_tsv("../annotationData/gencode.v24.annotation.hg19.middleTES.bed"          , col_names = FALSE)

colnames(tss) <- c("chr", "start", "end", "name", "score", "strand", "type", "symbol")
colnames(tes) <- c("chr", "start", "end", "name", "score", "strand", "type", "symbol")
tss$source <- "tss"
tes$source <- "tes"

genes <- rbind_list(tss, tes)

genes <- mutate(genes, width = genes$end - genes$start)

types <- rep("other", nrow(genes))
types[which(genes$type == "protein_coding")] <- "protein_coding"
types[which(genes$type == "processed_pseudogene")] <- "processed_pseudogene"
genes$type <- types

p <- ggplot(genes, aes(x = type, y = width, fill = source)) + geom_boxplot() + scale_y_log10() + theme_bw()
ggsave("plots/gene_type_width.pdf", p)

