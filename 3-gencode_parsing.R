# 2018-12-07 ----------
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz
gunzip gencode.v29.annotation.gtf.gz

gunzip gencode.v29.annotation.gff3.gz
liftOver -gff gencode.v29.annotation.gff3 ~/work/software/hg38ToHg19.over.chain.gz gencode.v29.annotation.hg19.gff3 gencode.v29.annotation.hg19.lostInTranslation.gff3

wc -l gencode.v29.*
2739404 gencode.v29.annotation.gff3
2735154 gencode.v29.annotation.hg19.gff3
8500 gencode.v29.annotation.hg19.lostInTranslation.gff3
252369 gencode.v29.transcripts.fa.gz

###
# R ------------
###
setwd("~/work/projects/cascade/data/annotation")
library(rtracklayer)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(furrr); plan(multisession(workers = 14))

t0 <- Sys.time()
gencode <- import("gencode.v29.annotation.hg19.gff3", format = "GFF", genome = "hg19")
Sys.time() - t0 # 1min

gencodeT <- as_tibble(as.data.frame(gencode, stringsAsFactors = FALSE))
transcript <- filter(gencodeT, type == "transcript")

getMiddleLineFor <- function(gene) {
    tempt <- dplyr::filter(transcript, gene_id == gene)
        if(tempt$strand[1] == "+") {
        tempt <- dplyr::arrange(tempt, start)
    } else if(tempt$strand[1] == "-") {
        tempt <- dplyr::arrange(tempt, end)
    }
    dplyr::slice(tempt, ceiling(nrow(tempt)/2))
}

t0 <- Sys.time()
midTranscript <- future_map_dfr(
    unique(transcript$gene_id),
    getMiddleLineFor
)
Sys.time() - t0 # 6 minutes

getFirstLineFor <- function(gene) {
    tempt <- dplyr::filter(transcript, gene_id == gene)
    if(tempt$strand[1] == "+") {
        tempt <- dplyr::arrange(tempt, start)
    } else if(tempt$strand[1] == "-") {
        tempt <- dplyr::arrange(tempt, end)
    }
    dplyr::slice(tempt, 1)
}

t0 <- Sys.time()
firstTranscript <- future_map_dfr(
    unique(transcript$gene_id),
    getFirstLineFor
)
Sys.time() - t0 # 6 minutes


midTranscript <- select(midTranscript, seqnames, start, end, gene_id,level, strand, gene_type, gene_name)
write.table(midTranscript, file = "gencode.v29.annotation.hg19.middleTSStranscript.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

midTSS <- midTranscript
midTSS$start <- ifelse(midTSS$strand == "+", midTSS$start, midTSS$end - 1)
midTSS$end <- midTSS$start + 1
write.table(midTSS, file = "gencode.v29.annotation.hg19.middleTSS.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(
    select(midTranscript, seqnames, start, end, gene_id,level, strand),
    file = "gencode.v29.annotation.hg19.middleTSStranscript.light.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    select(midTSS, seqnames, start, end, gene_id,level, strand),
    file = "gencode.v29.annotation.hg19.middleTSS.light.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

library(readr)
transcript <- read_tsv("gencode.v29.annotation.hg19.middleTSStranscript.light.bed", col_names = FALSE)
tss <- read_tsv("gencode.v29.annotation.hg19.middleTSS.light.bed", col_names = FALSE)

myChr <- paste0("chr", 1:22)
transcript <- filter(transcript, X1 %in% myChr)
tss <- filter(tss, X1 %in% myChr)

write.table(
    transcript,
    file = "gencode.v29.annotation.hg19.middleTSStranscript.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    tss,
    file = "gencode.v29.annotation.hg19.middleTSS.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

firstTranscript <- select(firstTranscript, seqnames, start, end, gene_id,level, strand, gene_type, gene_name)
firstTSS <- firstTranscript
firstTSS$start <- ifelse(firstTSS$strand == "+", firstTSS$start, firstTSS$end - 1)
firstTSS$end <- firstTSS$start + 1

write.table(
    filter(firstTranscript, seqnames %in% myChr),
    file = "gencode.v29.annotation.hg19.firstTSStranscript.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    filter(firstTSS, seqnames %in% myChr),
    file = "gencode.v29.annotation.hg19.firstTSS.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")



# TES ----------
getMiddleTESFor <- function(gene, reft) {
    tempt <- dplyr::filter(transcript, gene_id == gene)
    if(tempt$strand[1] == "+") {
        tempt <- dplyr::arrange(tempt, end)
    } else if(tempt$strand[1] == "-") {
        tempt <- dplyr::arrange(tempt, start)
    }
    dplyr::slice(tempt, ceiling(nrow(tempt)/2))
}


t0 <- Sys.time()
midTranscriptTES <- future_map_dfr(
    unique(transcript$gene_id),
    getMiddleTESFor
)
Sys.time() - t0 # 6 minutes


# revStrand <- ifelse(midTranscriptTES$strand == "+", "-", "+")
# midTranscriptTES$strand <- revStrand

midTranscriptTES <- select(midTranscriptTES, seqnames, start, end, gene_id,level, strand, gene_type, gene_name)

write.table(midTranscriptTES, file = "gencode.v29.annotation.hg19.middleTES.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    select(midTranscriptTES, seqnames, start, end, gene_id,level, strand),
    file = "gencode.v29.annotation.hg19.middleTES.light.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
myChr <- paste0("chr", 1:22)
midTranscriptTES <- filter(select(midTranscriptTES, seqnames, start, end, gene_id,level, strand), seqnames %in% myChr)
write.table(
    midTranscriptTES,
    file = "gencode.v29.annotation.hg19.middleTES.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

