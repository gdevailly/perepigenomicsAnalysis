setwd("/groups2/joshi_grp/guillaume/cascade/data/")

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(Repitools)
library(GenomicRanges)

source("../Rscripts/6-plotingFunctions.R")


# metadata loading, messy ---------------------------------------
metadata <- read_tsv("wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

myProms <- "../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.light.autosomes.bed"

roadmapExp <- list(
    pc = read_tsv("../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.pc"),
    nc = read_tsv("../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.nc"),
    rb = read_tsv("../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.rb")
)
roadmapExp <- do.call(rbind, roadmapExp)

refTable <- read_tsv(
    "../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.bed",
    col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)

metadata$id[!metadata$id %in% colnames(roadmapExp)] # missing E008, E017, E021, E022, check newer data?
metadata <- dplyr::filter(metadata, id %in% colnames(roadmapExp))

# build histone marks table --------------
hisFiles <- data_frame(
    file = list.files("histone/") %>%
        grep(".gz$", ., value = TRUE),
    name = gsub(".tagAlign.gz", "", file, fixed = TRUE)
)
hisFiles <- bind_cols(
    hisFiles,
    strsplit(hisFiles$name, "-") %>%
        map_df(~data_frame(cellCode = .x[1], ChIP = .x[2]))
)

# explore histone conditions -----------
table(hisFiles$ChIP)
table(hisFiles$cellCode, hisFiles$ChIP) # dark magic

filePath <- hisFiles$file[17]

# function definition ----------------------
importReadAlignAsGranges <- function(
    filePath,
    chromosomes_to_keep = paste0("chr", 1:22),
    col_types = "ciicic"
) {
    read_tsv(paste0("histone/", filePath), col_names = FALSE, col_types = col_types, progress = FALSE) %>% # import is long, file is huge
        dplyr::select(-X4, -X5) %>%
        dplyr::rename(chr = X1, start = X2, end = X3, strand = X6) %>%
        dplyr::filter(chr %in% chromosomes_to_keep) %>%
        with(., GRanges(chr, IRanges(start, end), strand = strand)) %>%
        return()
}

extractFeatureScoreMatrix <- function(
    genoRange,
    annoTable,
    up = 2500,
    down = 2500,
    freq = 100,
    ...
) {
    fs <- featureScores(
        genoRange,
        annoTable,
        up = up,
        down = down,
        freq = freq,
        use.strand = FALSE, # reads must not be on same strand as feature
        ...
    ) %>%
        tables %>%
        .[[1]]
    fs <- fs * 1E6 * 1000/freq # transform score into fpkm values
    return(fs)
}

# tests
t0 <- Sys.time() # 1 min 30
myGR <- importReadAlignAsGranges("E024-H3K4me3.tagAlign.gz")
Sys.time() - t0

# promoteur loading -----------------
tss_table <- read_tsv(myProms, col_names = FALSE, progress = FALSE)
colnames(tss_table) <- c("chr", "start", "end", "name", "score", "strand")
table(tss_table$strand)

t0 <- Sys.time()
myMat <- extractFeatureScoreMatrix(myGR, tss_table)
Sys.time() - t0

plotTraceAndKeyAndProfile(myMat) # \o/ works in pdf

myDF <- data.frame(
    gene_id = rownames(myMat),
    myMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))

myDF2 <- inner_join(
    select(roadmapExp, gene_id, E003),
    myDF,
    by = "gene_id"
) %>% arrange(desc(E003))

myMat2 <- select(myDF2, -gene_id, -E003) %>%
    as.matrix

pdf("test.pdf", width = 4, height = 10)
layout(matrix(1:3), heights = c(0.5, 0.15, 0.35))
plotTraceAndKeyAndProfile(myMat2)
dev.off()
