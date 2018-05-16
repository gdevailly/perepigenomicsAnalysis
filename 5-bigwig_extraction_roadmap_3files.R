# requires betdools:
# module load apps/gcc/BEDTools/2.25.0

setwd("/groups2/joshi_grp/guillaume/cascade/data/wgbs/roadmap/")
library(rtracklayer)
library(dplyr)
library(readr)
library(BSgenome.Hsapiens.UCSC.hg19)

# function -----------------
makeThreeBwFrom <- function(
    roadmapCode
) {
    message(paste("working on", roadmapCode))

    bw2bg <- paste0(
        "../../../../programmes/bigWigToBedGraph ",
        roadmapCode,
        "_WGBS_FractionalMethylation.bigwig ",
        roadmapCode,
        "_WGBS_FractionalMethylation.bedgraph"
    )
    message("bigwig to bedgraph...")
    system(bw2bg)
    message("done!")

    bedtools_map <- paste0(
        "bedtools map -a hg19_binned_w250_s100_s.bed -b ",
        roadmapCode,
        "_WGBS_FractionalMethylation.bedgraph -c 4,4,4 -o mean,count,sum -null NA > ",
        roadmapCode,
        "_WGBS_mean_count_sum.bed"
    )
    message("bedtools map...")
    system(bedtools_map)
    message("done!")

    message("processing bed file...")
    myFile <- read_tsv(paste0(roadmapCode ,"_WGBS_mean_count_sum.bed"), col_names = FALSE)
    colnames(myFile) <- c("chr", "start", "end", "mCpG_ratio", "CpG_sites", "mCpG_density")
    myFile <- filter(myFile, (!is.na(mCpG_ratio)) | (CpG_sites != 0) | (!is.na(mCpG_density)))

    myFile <- mutate(myFile, newStart = start + 75, newEnd = end - 76)
    myFile <- filter(myFile, newEnd > newStart) # we remove few windows at the end of each chr
    myFile_GRl <- list(
        mCpG_ratio   = with(myFile, GRanges(chr, IRanges(newStart, newEnd), score = mCpG_ratio  )),
        CpG_sites    = with(myFile, GRanges(chr, IRanges(newStart, newEnd), score = CpG_sites   )),
        mCpG_density = with(myFile, GRanges(chr, IRanges(newStart, newEnd), score = mCpG_density))
    )

    mySeqLen <- seqlengths(Hsapiens)
    mySeqLen <- mySeqLen[which(names(mySeqLen) %in% paste0("chr", 1:22))]
    mySeqLen <- mySeqLen[names(seqlengths(myFile_GRl[[1]]))]

    for (i in seq_along(myFile_GRl)) {
        seqlengths(myFile_GRl[[i]]) <- mySeqLen
    }

    for (i in seq_along(myFile_GRl)) {
        message(paste("Writing ", names(myFile_GRl)[i]))
        export(myFile_GRl[[i]], paste0(roadmapCode, "_", names(myFile_GRl)[i], "_hg19.bw"), format = "bw")
    }

    rm <- paste0(
        "rm ",
        roadmapCode,
        c("_WGBS_FractionalMethylation.bedgraph", "_WGBS_mean_count_sum.bed")
    )
    message("Cleaning files...")
    system(rm[1])
    system(rm[2])
    message(paste(roadmapCode, "done!"))
}

# code ----------------------

t0 <- Sys.time()
makeThreeBwFrom("E005")
Sys.time() - t0 # 10 min, works in screen, carefull with disk usage

files <- list.files()
myRMcodes <- grep("*_WGBS_FractionalMethylation.bigwig", files, value = TRUE) %>%
    gsub("_WGBS_FractionalMethylation.bigwig", "", ., fixed = TRUE)

t0 <- Sys.time()
for (i in 4:length(myRMcodes)) {
    makeThreeBwFrom(myRMcodes[i])
    Sys.time() - t0
}

# expanding coverage bigwig for concistency ----------------
makeSootherBigWigCoverage <- function(
    roadmapCode
) {
    message(paste("working on", roadmapCode))

    bw2bg <- paste0(
        "../../../../programmes/bigWigToBedGraph ",
        roadmapCode,
        "_WGBS_ReadCoverage.bigwig ",
        roadmapCode,
        "_WGBS_ReadCoverage.bedgraph"
    )
    message("bigwig to bedgraph...")
    system(bw2bg)
    message("done!")

    bedtools_map <- paste0(
        "bedtools map -a hg19_binned_w250_s100_s.bed -b ",
        roadmapCode,
        "_WGBS_ReadCoverage.bedgraph -c 4 -o mean -null NA > ",
        roadmapCode,
        "_WGBS_ReadCoverage_mean.bed"
    )
    message("bedtools map...")
    system(bedtools_map)
    message("done!")

    message("processing bed file...")
    myFile <- read_tsv(paste0(roadmapCode ,"_WGBS_ReadCoverage_mean.bed"), col_names = FALSE)
    colnames(myFile) <- c("chr", "start", "end", "coverage")
    myFile <- filter(myFile, (!is.na(coverage)))

    myFile <- mutate(myFile, newStart = start + 75, newEnd = end - 76) # centering the window
    myFile <- filter(myFile, newEnd > newStart) # we remove few windows at the end of each chr
    myFile_GRl <- with(myFile, GRanges(chr, IRanges(newStart, newEnd), score = as.numeric(coverage)))

    mySeqLen <- seqlengths(Hsapiens)
    mySeqLen <- mySeqLen[which(names(mySeqLen) %in% paste0("chr", 1:22))]
    mySeqLen <- mySeqLen[names(seqlengths(myFile_GRl))]

    seqlengths(myFile_GRl) <- mySeqLen

    message(paste("Writing bigwig coverage file"))
    export(myFile_GRl, paste0(roadmapCode, "_coverage_hg19.bw"), format = "bw")

    rm <- paste0(
        "rm ",
        roadmapCode,
        c("_WGBS_ReadCoverage.bedgraph", "_WGBS_ReadCoverage_mean.bed")
    )
    message("Cleaning files...")
    system(rm[1])
    system(rm[2])
    message(paste(roadmapCode, "done!"))
}

files <- list.files()
myRMcodes <- grep("*_WGBS_FractionalMethylation.bigwig", files, value = TRUE) %>%
    gsub("_WGBS_FractionalMethylation.bigwig", "", ., fixed = TRUE)

t0 <- Sys.time()
for (i in 2:length(myRMcodes)) {
    makeSootherBigWigCoverage(myRMcodes[i])
    Sys.time() - t0
}



