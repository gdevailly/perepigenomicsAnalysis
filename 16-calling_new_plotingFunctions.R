setwd("/groups2/joshi_grp/guillaume/cascade/data/wgbs/roadmap")
library(dplyr)
library(plotrix)
library(seqplots)
library(readr)
library(svglite)

source("../../../Rscripts/6-plotingFunctions.R")

metadata <- read_tsv("EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

myProms <- "../../../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.light.autosomes.bed"

roadmapExp <- list(
    pc = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.pc"),
    nc = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.nc"),
    rb = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.rb")
)
roadmapExp <- do.call(rbind, roadmapExp)

refTable <- read_tsv(
    "../../../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.bed",
    col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)

metadata$id[!metadata$id %in% colnames(roadmapExp)] # missing E008, E017, E021, E022, check newer data?
metadata <- filter(metadata, id %in% colnames(roadmapExp))

# TSS ------------
t0 <- Sys.time() # 45 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataFor(
        metadata$id[i],
        myProms,
        roadmapExp,
        refgenome = "hg19",
        bin = 100L,
        rm0 = TRUE,
        xmin = 2500L, xmax = 2500L, type = "pf",
        add_heatmap = TRUE,
        verbose = TRUE
    )
    dataForPlot <- addGeneTypeInfo(dataForPlot, refTable)
    geneType <- rep("other", nrow(dataForPlot))
    geneType[which(dataForPlot$gene_type == "protein_coding")] <- "protein_coding"
    geneType[which(dataForPlot$gene_type == "processed_pseudogene")] <- "processed_pseudogene"
    dataForPlot2gt <- dataForPlot
    dataForPlot2gt$gene_type <- geneType
    dataForPlotPC <- filter(dataForPlot, gene_type == "protein_coding")
    dataForPlotPG <- filter(dataForPlot, gene_type == "processed_pseudogene")
    dataForPlotOther <- filter(dataForPlot, !gene_type %in% c("protein_coding", "processed_pseudogene"))
    geneType <- rep("other", nrow(dataForPlotOther))
    geneType[which(dataForPlotOther$gene_type == "lincRNA")] <- "lincRNA"
    geneType[which(dataForPlotOther$gene_type == "antisense")] <- "antisense"
    geneType[which(dataForPlotOther$gene_type == "miRNA")] <- "miRNA"
    dataForPlotOther$gene_type <- geneType

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
        file = paste0("../../../plots/profileHM/tss/", metadata$id[i], "_", metadata$short[i], "_tss_profile.pdf"),
        width = 8.3, height = 5.8
    )
    plotMetricAndProfile(
        dataForPlot2gt,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("grey", "red", "yellow")),
        naColour = "yellow",
        tints  = c("red", "blue", "purple", "sienna4")
    )
    plotMetricAndProfile(
        dataForPlotPC,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("yellow")),
        naColour = "yellow",
        tints  = c("red", "blue", "purple", "sienna4")
    )
    plotMetricAndProfile(
        dataForPlotPG,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("red")),
        naColour = "yellow",
        tints  = c("red", "blue", "purple", "sienna4")
    )
    plotMetricAndProfile(
        dataForPlotOther,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        naColour = "yellow",
        tints  = c("red", "blue", "purple", "sienna4")
    )
    dev.off()
    message(paste(metadata$id[i], "done!"))
    message(Sys.time() - t0)
}


# TES ----------------------
TES <- "../../../../annotationData/gencode.v24.annotation.hg19.middleTES.light.autosomes.bed"

t0 <- Sys.time() # 45 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataFor(
        metadata$id[i],
        TES,
        roadmapExp,
        refgenome = "hg19",
        bin = 100L,
        rm0 = TRUE,
        xmin = 2500L, xmax = 2500L, type = "ef",
        add_heatmap = TRUE,
        verbose = TRUE
    )
    dataForPlot <- addGeneTypeInfo(dataForPlot, refTable)
    geneType <- rep("other", nrow(dataForPlot))
    geneType[which(dataForPlot$gene_type == "protein_coding")] <- "protein_coding"
    geneType[which(dataForPlot$gene_type == "processed_pseudogene")] <- "processed_pseudogene"
    dataForPlot2gt <- dataForPlot
    dataForPlot2gt$gene_type <- geneType
    dataForPlotPC <- filter(dataForPlot, gene_type == "protein_coding")
    dataForPlotPG <- filter(dataForPlot, gene_type == "processed_pseudogene")
    dataForPlotOther <- filter(dataForPlot, !gene_type %in% c("protein_coding", "processed_pseudogene"))
    geneType <- rep("other", nrow(dataForPlotOther))
    geneType[which(dataForPlotOther$gene_type == "lincRNA")] <- "lincRNA"
    geneType[which(dataForPlotOther$gene_type == "antisense")] <- "antisense"
    geneType[which(dataForPlotOther$gene_type == "miRNA")] <- "miRNA"
    dataForPlotOther$gene_type <- geneType

    zlim1 <- c(0, 8)
    zlim4 <- c(10, 50)
    if(metadata$id[i] %in% c(
        "E005", "E007", "E011", "E012", "E013", "E016", "E065", "E066", "E070", "E071", "E095", "E106"
    )) zlim4 <- c(30, 100)
    if(metadata$id[i] %in% c("E084", "E085")) {
        zlim1 <- c(0, 4)
        zlim4 <- c(0, 20)
    }
    zlims <- list(zlim1, c(0,1), c(1, 4), zlim4)

    pdf(
        file = paste0("../../../plots/profileHM/tes/", metadata$id[i], "_", metadata$short[i], "_tes_profile.pdf"),
        width = 8.3, height = 5.8
    )
    plotMetricAndProfile(
        dataForPlot2gt,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("grey", "red", "yellow")),
        naColour = "yellow",
        xlabels = c("-2.5kb", "TES", "+2.5kb"),
        tints  = c("red", "blue", "purple", "sienna4")
    )
    plotMetricAndProfile(
        dataForPlotPC,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("yellow")),
        naColour = "yellow",
        xlabels = c("-2.5kb", "TES", "+2.5kb"),
        tints  = c("red", "blue", "purple", "sienna4")
    )
    plotMetricAndProfile(
        dataForPlotPG,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("red")),
        naColour = "yellow",
        xlabels = c("-2.5kb", "TES", "+2.5kb"),
        tints  = c("red", "blue", "purple", "sienna4")
    )
    plotMetricAndProfile(
        dataForPlotOther,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        naColour = "yellow",
        xlabels = c("-2.5kb", "TES", "+2.5kb"),
        tints  = c("red", "blue", "purple", "sienna4")
    )
    dev.off()
    message(paste(metadata$id[i], "done!"))
    message(Sys.time() - t0)
}


