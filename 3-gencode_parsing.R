# 2016-05-06 ----------
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz
gunzip gencode.v24.annotation.gtf.gz
../../liftOver -gff gencode.v24.annotation.gtf hg38ToHg19.over.chain.gz gencode.v24.annotation.hg19.gtf gencode.v24.annotation.hg19.lostInTranslation.gtf
wc -l gencode.v24.*
2572845 gencode.v24.annotation.gtf
2568520 gencode.v24.annotation.hg19.gtf
8650 gencode.v24.annotation.hg19.lostInTranslation.gtf

###
# R ------------
###
setwd("/groups2/joshi_grp/guillaume/annotationData")
library(rtracklayer)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(parallel)

t0 <- Sys.time()
gencode <- import("gencode.v24.annotation.hg19.gtf", format = "gtf", genome = "hg19")
Sys.time() - t0 # 3min

gencodeT <- tbl_df(as.data.frame(gencode, stringsAsFactors = FALSE))
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

midTranscript <- do.call(rbind, mclapply(
    unique(transcript$gene_id),
    getMiddleLineFor,
    mc.cores = 24
))
Sys.time() - t0 # mostly slowish rbind. dplyr solution? 25 min


midTranscript <- select(midTranscript, seqnames, start, end, gene_id,level, strand, gene_type, gene_name)
write.table(midTranscript, file = "gencode.v24.annotation.hg19.middleTSStranscript.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

midTSS <- midTranscript
midTSS$start <- ifelse(midTSS$strand == "+", midTSS$start, midTSS$end - 1)
midTSS$end <- midTSS$start + 1
write.table(midTSS, file = "gencode.v24.annotation.hg19.middleTSS.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(
    select(midTranscript, seqnames, start, end, gene_id,level, strand),
    file = "gencode.v24.annotation.hg19.middleTSStranscript.light.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    select(midTSS, seqnames, start, end, gene_id,level, strand),
    file = "gencode.v24.annotation.hg19.middleTSS.light.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

library(readr)
transcript <- read_tsv("gencode.v24.annotation.hg19.middleTSStranscript.light.bed", col_names = FALSE)
tss <- read_tsv("gencode.v24.annotation.hg19.middleTSS.light.bed", col_names = FALSE)

myChr <- paste0("chr", 1:22)
transcript <- filter(transcript, X1 %in% myChr)
tss <- filter(tss, X1 %in% myChr)

write.table(
    transcript,
    file = "gencode.v24.annotation.hg19.middleTSStranscript.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    tss,
    file = "gencode.v24.annotation.hg19.middleTSS.light.autosomes.bed",
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
midTranscriptTES <- do.call(rbind, mclapply(
    unique(transcript$gene_id),
    getMiddleTESFor,
    mc.cores = 6
))
Sys.time() - t0 # 33 min mostly slowish rbind.

# revStrand <- ifelse(midTranscriptTES$strand == "+", "-", "+")
# midTranscriptTES$strand <- revStrand

midTranscriptTES <- select(midTranscriptTES, seqnames, start, end, gene_id,level, strand, gene_type, gene_name)

write.table(midTranscriptTES, file = "gencode.v24.annotation.hg19.middleTES.bed",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(
    select(midTranscriptTES, seqnames, start, end, gene_id,level, strand),
    file = "gencode.v24.annotation.hg19.middleTES.light.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
myChr <- paste0("chr", 1:22)
midTranscriptTES <- filter(select(midTranscriptTES, seqnames, start, end, gene_id,level, strand), seqnames %in% myChr)
write.table(
    midTranscriptTES,
    file = "gencode.v24.annotation.hg19.middleTES.light.autosomes.bed",
    quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

