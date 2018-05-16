setwd("/groups2/joshi_grp/guillaume/cascade/data/")
source("../Rscripts/6-plotingFunctions.R")
source("../Rscripts/20-functions_for_histoneMarks.R")

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(parallel)

metadata <- read_tsv("wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

load("Rdata/innerExonQuant.RData")
load("Rdata/innerExonRatio.RData")

metadata$id[!metadata$id %in% colnames(innerExonRatio)] # missing E008, E017, E021, E022
metadata <- filter(metadata, id %in% colnames(innerExonRatio))

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

his_md <- left_join(
    dplyr::filter(hisFiles, !ChIP %in% c("DNase", "DNaseControl", "Input")),
    dplyr::filter(hisFiles, ChIP == "Input") %>%
        dplyr::select(file, cellCode) %>%
        dplyr::rename(input_file = file),
    by = "cellCode"
)

# DNAseI -----------------------------
DNAse_md <- inner_join(
    dplyr::filter(hisFiles, ChIP  == "DNase") %>%
        dplyr::select(-name, -ChIP),
    dplyr::filter(hisFiles, ChIP  == "DNaseControl") %>%
        dplyr::select(-name, -ChIP),
    by = "cellCode"
) %>% dplyr::rename(DNAse_file = file.x, Control_file = file.y)


# exon FPKM plot ------------------
exonInfoPath <- "inner_exon.bed"
exonTable <- read_tsv(exonInfoPath, col_names = FALSE, progress = FALSE)
colnames(exonTable) <- c("chr", "start", "end", "name", "score", "strand")
innerExonQuant <- dplyr::select(innerExonQuant, -name) %>%
    dplyr::rename(name = exon_location)

# png versions -------------------

plotHisModDataForLine <- function(i) {
    hisModData <- prepareHisModDataExons(
        i, annoTable = exonTable, metadataTable = his_md, metricTable = innerExonQuant,
        up = 1000, down = 1000, freq = 50
    )

    hisModName <- his_md[i, ]$ChIP
    cellID <-  his_md[i, ]$cellCode

    plotMetricAndProfile(
        hisModData,
        prefix = c("HisMod"  , "Input"),
        labels = c(hisModName, "Input"),
        tints  = c("sienna1" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        withGeneType = FALSE,
        xlabels = c("-1kb", "Exon start", "+1kb"),
        npix_height = 600
    )
    message(paste(his_md[i, ]$name, "is done!"))
}

t0 <- Sys.time() # 6h20
mclapply(
    seq_len(nrow(his_md)),
    function(i) {
        hisMark <- his_md[i, ]$ChIP
        cellID  <- his_md[i, ]$cellCode
        if(!file.exists(paste0("../appPlots/exons/", hisMark))){
            dir.create(paste0("../appPlots/exons/", hisMark))
            dir.create(paste0("../appPlots/exons/", hisMark, "/byFpkm"))
        }
        png(
            file = paste0(
                "../appPlots/exons/", hisMark, "/byFpkm/", cellID, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotHisModDataForLine(i)
        dev.off()
    },
    mc.cores = 8 # big memmory footprint
) %>% invisible()
Sys.time() - t0

# DNAseI -----------------------
plotDNAseDataForLine <- function(i) {
    dnaseData <- prepareDNAseDataExons(
        DNAse_md$cellCode[i], annoTable = exonTable, metadataTable = DNAse_md, metricTable = innerExonQuant,
        up = 1000, down = 1000, freq = 50
    )

    cellID <- DNAse_md$cellCode[i]

    plotMetricAndProfile(
        dnaseData,
        prefix = c("DNAse"  , "Control"),
        labels = c("DNAse"  , "Control"),
        tints  = c("dodgerblue" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        withGeneType = FALSE,
        xlabels = c("-1kb", "Exon start", "+1kb"),
        npix_height = 600
    )
    message(paste(cellID, "is done!"))
}

t0 <- Sys.time() # 31 minutes
mclapply(
    seq_len(nrow(DNAse_md)),
    function(i) {
        cellID <- DNAse_md$cellCode[i]
        if(!file.exists(paste0("../appPlots/exons/dnase1"))){
            dir.create(paste0("../appPlots/exons/dnase1"))
            dir.create(paste0("../appPlots/exons/dnase1/byFpkm"))
        }
        png(
            file = paste0(
                "../appPlots/exons/dnase1/byFpkm/", cellID, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotDNAseDataForLine(i)
        dev.off()
    },
    mc.cores = 7 # big memmory footprint
) %>% invisible()
Sys.time() - t0


# exon ratio ----------------------

innerExonRatio <- select(innerExonRatio, -name) %>%
    dplyr::rename(name = exon_location)

plotHisModDataForLine <- function(i) {

    hisModData <- prepareHisModDataExons(
        i, annoTable = exonTable, metadataTable = his_md, metricTable = innerExonRatio,
        up = 1000, down = 1000, freq = 50
    )

    hisModName <- his_md[i, ]$ChIP
    cellID <-  his_md[i, ]$cellCode

    bins <- nrow(hisModData)
    bins[which(hisModData$exp >= 3)] <- 1
    bins[which(hisModData$exp < 3 & hisModData$exp >= 1)] <- 2
    bins[which(hisModData$exp < 1 & hisModData$exp >= 1/3)] <- 3
    bins[which(hisModData$exp < 1/3)] <- 4

    plotMetricAndProfile(
        hisModData,
        prefix = c("HisMod"  , "Input"),
        labels = c(hisModName, "Input"),
        tints  = c("blue" , "grey"),
        zlims  = list(c(0, 1) , c(0, 1)),
        title = paste(hisModName, gsub("_", " ", filter(metadata, id == cellID)$name, fixed = TRUE)),
        withGeneType = FALSE,
        xlabels = c("-1kb", "Exon start", "+1kb"),
        expLim = c(-4, 4),
        expTransf = function(x) log2(x+0.01),
        bins = bins,
        exp.text = "log2(exon/gene)"
    )

    message(paste(his_md[i, ]$name, "is done!"))
}

t0 <- Sys.time() # 11.62 hours with 4 cores
mclapply(
    seq_len(nrow(his_md)),
    function(i) {
        hisMark <- his_md[i, ]$ChIP
        cellID <-  his_md[i, ]$cellCode
        if(!file.exists(paste0("../plots/profileHM/exons/byRatio/", hisMark))){
            dir.create(paste0("../plots/profileHM/exons/byRatio/", hisMark))
        }
        pdf(
            file = paste0(
                "../plots/profileHM/exons/byRatio/", hisMark, "/",  hisMark, "_",
                dplyr::filter(metadata, id == cellID)$id, "_", dplyr::filter(metadata, id == cellID)$short,
                "_exons_fpkm_profile.pdf"
            ),
            width = 5.3, height = 5.8
        )
        plotHisModDataForLine(i) %>% invisible
        dev.off()
    },
    mc.cores = 4 # big memmory footprint
) %>% invisible
Sys.time() - t0

# exon width ---------------------

plotHisModDataForLine <- function(i) {

    hisModData <- prepareHisModDataExons(
        i, annoTable = exonTable, metadataTable = his_md, metricTable = innerExonRatio,
        up = 1000, down = 1000, freq = 50
    )

    hisModData <- mutate(hisModData, exp = gsub("chr[[:digit:]]+:", "", name) %>%
        gsub("<[-]?1$", "", .)  %>%
        strsplit("-") %>%
        map_int(~as.integer(dplyr::last(.x)) - as.integer(dplyr::first(.x)))
    ) %>%
        arrange(desc(exp))

    hisModName <- his_md[i, ]$ChIP
    cellID <-  his_md[i, ]$cellCode

    plotMetricAndProfile(
        hisModData,
        prefix = c("HisMod"  , "Input"),
        labels = c(hisModName, "Input"),
        tints  = c("blue" , "grey"),
        zlims  = list(c(0, 1) , c(0, 1)),
        title = paste(hisModName, gsub("_", " ", filter(metadata, id == cellID)$name, fixed = TRUE)),
        withGeneType = FALSE,
        xlabels = c("-1kb", "Exon start", "+1kb"),
        expLim = c(0, 1000),
        expTransf = function(x) x,
        exp.text = "exon width"
    )

    message(paste(his_md[i, ]$name, "is done!"))
}

t0 <- Sys.time() # 11 hours with 4 cores
mclapply(
    seq_len(nrow(his_md)),
    function(i) {
        hisMark <- his_md[i, ]$ChIP
        cellID <-  his_md[i, ]$cellCode
        if(!file.exists(paste0("../plots/profileHM/exons/byWidth/", hisMark))){
            dir.create(paste0("../plots/profileHM/exons/byWidth/", hisMark))
        }
        pdf(
            file = paste0(
                "../plots/profileHM/exons/byWidth/", hisMark, "/",  hisMark, "_",
                dplyr::filter(metadata, id == cellID)$id, "_", dplyr::filter(metadata, id == cellID)$short,
                "_exons_fpkm_profile.pdf"
            ),
            width = 5.3, height = 5.8
        )
        plotHisModDataForLine(i) %>% invisible
        dev.off()
    },
    mc.cores = 4 # big memmory footprint
) %>% invisible
Sys.time() - t0


