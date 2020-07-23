library(here)
library(rtracklayer)
library(tidyverse)
library(parallel)

t0 <- Sys.time()
gencode <- import(here("data", "annotation", "gencode.v24.annotation.hg19.gtf"), format = "gtf")
Sys.time() - t0 # 1min

gencodeT <- as_tibble(gencode)

unique(gencodeT$type)

gene       <- filter(gencodeT, type == "gene"       & gene_type == "protein_coding")
transcript <- filter(gencodeT, type == "transcript" & gene_type == "protein_coding")
exon       <- filter(gencodeT, type == "exon"       & gene_type == "protein_coding")

unique(exon$gene_type)

# take a gene name (ensembl), a transcript table and exon table.
# return middle exons (not first, not start), of the shorter transcripts
getMiddleExonsFor <- function(myGene, transcript, exon) {
    myTranscripts <- filter(transcript, gene_id == myGene) %>% arrange(width)
    myExons <- filter(exon, transcript_id == myTranscripts$transcript_id[1])
    # trust gencode sorting, may need some additional sorting?
    return(myExons[c(-1,-nrow(myExons)), ])
}

getMiddleExonsFor("ENSG00000049249.8", transcript, exon)

t0 <- Sys.time()
myExons <- bind_rows(
    mclapply(
        gene$gene_id,
        function(x) getMiddleExonsFor(x, transcript, exon),
        mc.cores = 10
    )
)
Sys.time() - t0 # 1.2 minutes, single threaded

library(ggplot2)
ggplot(myExons, aes(width)) + geom_density() + scale_x_log10() + theme_bw() + ggtitle("Exons length")
p <- ggplot(myExons, aes(width)) + geom_density() + coord_cartesian(x = c(0, 500)) + theme_bw() + ggtitle("Exons length")
setwd("/groups2/joshi_grp/guillaume/cascade/")
ggsave("plots/exons_width.pdf", p)

ggplot(myExons, aes(gene_id)) + geom_bar() + theme_bw()

p <- ggplot(data.frame(table(myExons$gene_id)), aes(Freq)) + geom_bar() +  coord_cartesian(x = c(0, 20)) + theme_bw() +
        labs(title = "Middle exons per gene", x = "Exons per gene")
ggsave("plots/exons_per_gene.pdf", p)


myExonsLight <- select(myExons, seqnames, start, end, gene_id, level, strand, gene_type, gene_name)

setwd("/groups2/joshi_grp/guillaume/annotationData")
write.table(myExonsLight, file = "gencode.v24.annotation.hg19.middle.exons.pc.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    select(myExonsLight, seqnames, start, end, gene_id,level, strand),
    file = "gencode.v24.annotation.hg19.middle.exons.pc.light.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
myChr <- paste0("chr", 1:22)
myExonsLightAuto <- filter(select(myExonsLight, seqnames, start, end, gene_id,level, strand), seqnames %in% myChr)
write.table(
    myExonsLightAuto,
    file = "gencode.v24.annotation.hg19.middle.exons.pc.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")





