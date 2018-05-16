setwd("/groups2/joshi_grp/guillaume/cascade/data/wgbs/roadmap")
library(dplyr)
library(readr)
library(ggplot2)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)

t0 <- Sys.time()
E003 <- read_tsv("E003_WGBS_mean_count_sum.bed", col_names = FALSE)
Sys.time() - t0 # 30 sec
colnames(E003) <- c("chr", "start", "end", "mCpG_ratio", "CpG_sites", "mCpG_density")

nrow(E003) # 28.810.343
E003 <- filter(E003, (!is.na(mCpG_ratio)) | (CpG_sites != 0) | (!is.na(mCpG_density)))
nrow(E003) # 20.385.296

E003 <- mutate(E003, newStart = start + 75, newEnd = end - 76)
E003 <- filter(E003, newEnd > newStart) # we remove few windows at the end of each chr
E003_GRl <- list(
    mCpG_ratio = with(E003, GRanges(chr, IRanges(newStart, newEnd), score = mCpG_ratio)),
    CpG_sites = with(E003, GRanges(chr, IRanges(newStart, newEnd), score = CpG_sites)),
    mCpG_density = with(E003, GRanges(chr, IRanges(newStart, newEnd), score = mCpG_density))
)

mySeqLen <- seqlengths(Hsapiens)
mySeqLen <- mySeqLen[which(names(mySeqLen) %in% paste0("chr", 1:22))]
mySeqLen <- mySeqLen[names(seqlengths(E003_GRl[[1]]))]

for (i in seq_along(E003_GRl)) {
    seqlengths(E003_GRl[[i]]) <- mySeqLen
}

t0 <- Sys.time()
for (i in seq_along(E003_GRl)) {
    export(E003_GRl[[i]], paste0("E003_", names(E003_GRl[i]), "_hg19.bw"), format = "bw")
}
Sys.time() - t0 # 2min



