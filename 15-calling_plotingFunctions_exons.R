setwd("/groups2/joshi_grp/guillaume/cascade/data/")
source("../Rscripts/6-plotingFunctions.R")

library(dplyr)
library(plotrix)
library(seqplots)
library(readr)
library(svglite)

metadata <- read_tsv("wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

load("Rdata/innerExonQuant.RData")
load("Rdata/innerExonRatio.RData")

metadata$id[!metadata$id %in% colnames(innerExonRatio)] # missing E008, E017, E021, E022
metadata <- filter(metadata, id %in% colnames(innerExonRatio))

exonInfo <- innerExonQuant[, c("chr", "start", "end", "exon_location", "score", "strand")]
write.table(exonInfo, file = "inner_exon.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

setwd("wgbs/roadmap")

# exon FPKM plot ------------------
exonInfoPath <- "../../inner_exon.bed"

t0 <- Sys.time() # 45 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataForExons(
        metadata$id[i],
        exonInfoPath,
        innerExonQuant,
        refgenome = "hg19",
        bin = 50L,
        rm0 = TRUE,
        xmin = 1000L, xmax = 1000L, type = "mf",
        add_heatmap = TRUE,
        verbose = TRUE
    )

    zlim1 <- c(0, 20)
    zlim4 <- c(10, 50)
    if(metadata$id[i] %in% c(
        "E005", "E007", "E011", "E012", "E013", "E016", "E065", "E066", "E070", "E071", "E095", "E106"
    )) zlim4 <- c(30, 100)
    if(metadata$id[i] %in% c("E084", "E085")) {
        zlim1 <- c(0, 10)
        zlim4 <- c(0, 20)
    }
    zlims <- list(zlim1, c(0,1), c(1, 4), zlim4)

    pdf(
        file = paste0("../../../plots/profileHM/exons/byFPKM/", metadata$id[i], "_", metadata$short[i], "_exonsByFPKM_profile.pdf"),
        width = 8.3, height = 5.8
    )
    plotMetricAndProfile(
        dataForPlot,
        nbin = 5,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("yellow")),
        naColour = "yellow",
        withGeneType = FALSE,
        expTransf = function(x) log10(x + 1),
        xlabels = c("-1kb", "center", "+1kb"),
        expLim = c(0, 5),
        exp.text = "log10(FPKM + 1)",
        reversedZOrder = TRUE,
        tints  = c("red"        , "blue"      , "purple"      , "sienna4"   )
    )
    dev.off()

    message(paste(metadata$id[i], "done!"))
    message(Sys.time() - t0)
}

# png version --------------------
t0 <- Sys.time() # 45 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataForExons(
        metadata$id[i],
        exonInfoPath,
        innerExonQuant,
        refgenome = "hg19",
        bin = 50L,
        rm0 = TRUE,
        xmin = 1000L, xmax = 1000L, type = "pf",
        add_heatmap = TRUE,
        verbose = TRUE
    )

    zlim1 <- c(0, 20)
    zlim4 <- c(10, 50)
    if(metadata$id[i] %in% c(
        "E005", "E007", "E011", "E012", "E013", "E016", "E065", "E066", "E070", "E071", "E095", "E106"
    )) zlim4 <- c(30, 100)
    if(metadata$id[i] %in% c("E084", "E085")) {
        zlim1 <- c(0, 10)
        zlim4 <- c(0, 20)
    }
    zlims <- list(zlim1, c(0,1), c(1, 4), zlim4)

    png(
        file = paste0("../../../appPlots/exons/wgbs/byFpkm/", metadata$id[i], ".png"),
        width = 8.3, height = 5.8, units = "in", res = 300, pointsize = 13
    )
    plotMetricAndProfile(
        dataForPlot,
        nbin = 5,
        zlims  = zlims,
        raster = TRUE,
        title = "",
        withGeneType = FALSE,
        expTransf = function(x) log10(x + 1),
        xlabels = c("-1kb", "Exon start", "+1kb"),
        expLim = c(0, 4),
        exp.text = "log10(FPKM + 1)",
        reversedZOrder = TRUE,
        tints  = c("red"        , "blue"      , "purple"      , "sienna4"   )
    )
    dev.off()

    message(paste(metadata$id[i], "done!"))
}
Sys.time() - t0

# exon ratio plot ------------------
exonInfoPath <- "../../inner_exon.bed"

t0 <- Sys.time() # 12 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataForExons(
        metadata$id[i],
        exonInfoPath,
        innerExonRatio,
        refgenome = "hg19",
        bin = 50L,
        rm0 = TRUE,
        xmin = 1000L, xmax = 1000L, type = "mf",
        add_heatmap = TRUE,
        verbose = TRUE
    )

    zlim1 <- c(0, 20)
    zlim4 <- c(10, 50)
    if(metadata$id[i] %in% c(
        "E005", "E007", "E011", "E012", "E013", "E016", "E065", "E066", "E070", "E071", "E095", "E106"
    )) zlim4 <- c(30, 100)
    if(metadata$id[i] %in% c("E084", "E085")) {
        zlim1 <- c(0, 10)
        zlim4 <- c(0, 20)
    }
    zlims <- list(zlim1, c(0,1), c(1, 4), zlim4)

    bins <- numeric(21125)
    bins[which(dataForPlot$exp >= 3)] <- 1
    bins[which(dataForPlot$exp < 3 & dataForPlot$exp >= 1)] <- 2
    bins[which(dataForPlot$exp < 1 & dataForPlot$exp >= 1/3)] <- 3
    bins[which(dataForPlot$exp < 1/3)] <- 4

    pdf(
        file = paste0("../../../plots/profileHM/exons/byRatio/", metadata$id[i], "_", metadata$short[i], "_exonsByRatio_profile_yellow_1.58_n.pdf"),
        width = 8.3, height = 5.8
    )
    plotMetricAndProfile(
        dataForPlot,
        bins = bins,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("yellow")),
        naColour = "yellow",
        withGeneType = FALSE,
        expTransf = function(x) log2(x+0.01),
        xlabels = c("-1kb", "center", "+1kb"),
        expLim = c(-5, 5),
        exp.text = "log2(exon/gene)",
        reversedZOrder = TRUE,
        abline.v = 0,
        tints  = c("red"        , "blue"      , "purple"      , "sienna4"   )
    )
    dev.off()

    message(paste(metadata$id[i], "done!"))
    message(Sys.time() - t0)

} # 12 minutes

## exon length vs exon ratio bins
library(cowplot)
exonInfoPath <- "../../inner_exon.bed"

t0 <- Sys.time() # 22 sec
pdf(
    file = "../../../plots/profileHM/exons/byRatio/exon_width_by_ratio_bins.pdf",
    width = 4, height = 4
)
for(i in seq_len(nrow(metadata))){
    ewt <- select_(innerExonRatio, "chr", "start", "end", "exon_location", exp = metadata$id[i]) %>%
        mutate(width = end - start)
    bins <- numeric(nrow(ewt))
    bins[which(ewt$exp >= 3)] <- 1
    bins[which(ewt$exp < 3 & ewt$exp >= 1)] <- 2
    bins[which(ewt$exp < 1 & ewt$exp >= 1/3)] <- 3
    bins[which(ewt$exp < 1/3)] <- 4
    ewt$bins <- factor(bins)

    p <- ggplot(ewt, aes(x = bins, y = width, fill = bins)) +
        geom_boxplot() +
        scale_y_log10() +
        labs(title = gsub("_", " ", metadata$name[i], fixed = TRUE)) +
        theme(legend.position = "none")
    print(p)
}
dev.off()
Sys.time() - t0

## ----------------
bins <- numeric(21125)
bins[which(dataForPlot$exp >= 2)] <- 1
bins[which(dataForPlot$exp < 2 & dataForPlot$exp >= 1)] <- 2
bins[which(dataForPlot$exp < 1 & dataForPlot$exp >= 0.5)] <- 3
bins[which(dataForPlot$exp < 0.5)] <- 4

pdf(
    file = paste0("../../../plots/profileHM/exons/byRatio/", metadata$id[i], "_", metadata$short[i], "_exonsByRatio_profile_bins_newer3.pdf"),
    width = 8.3, height = 5.8
)
plotMetricAndProfile(
    dataForPlot,
    bins = bins,
    zlims  = zlims,
    raster = TRUE,
    title = gsub("_", " ", metadata$name[i], fixed = TRUE),
    geneTypePalette = colorRampPalette(c("yellow")),
    naColour = "grey",
    withGeneType = FALSE,
    expTransf = function(x) log10(x+1),
    xlabels = c("-1kb", "center", "+1kb"),
    expLim = c(-5, 5),
    exp.text = "log10(exon FPKM + 1)",
    reversedZOrder = TRUE,
    abline.v = c(-1, 0, 1)
)
dev.off()

