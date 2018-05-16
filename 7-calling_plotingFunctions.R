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
    plotExpressionAndProfile(
        dataForPlot2gt,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("grey", "red", "yellow")),
        naColour = "white"
    )
    plotExpressionAndProfile(
        dataForPlotPC,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("yellow")),
        naColour = "white"
    )
    plotExpressionAndProfile(
        dataForPlotPG,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("red")),
        naColour = "white"
    )
    plotExpressionAndProfile(
        dataForPlotOther,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        naColour = "white"
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
    plotExpressionAndProfile(
        dataForPlot2gt,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("grey", "red", "yellow")),
        naColour = "white",
        xlabels = c("-2.5kb", "TES", "+2.5kb")
    )
    plotExpressionAndProfile(
        dataForPlotPC,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("yellow")),
        naColour = "white",
        xlabels = c("-2.5kb", "TES", "+2.5kb")
    )
    plotExpressionAndProfile(
        dataForPlotPG,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("red")),
        naColour = "white",
        xlabels = c("-2.5kb", "TES", "+2.5kb")
    )
    plotExpressionAndProfile(
        dataForPlotOther,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        naColour = "white",
        xlabels = c("-2.5kb", "TES", "+2.5kb")
    )
    dev.off()
    message(paste(metadata$id[i], "done!"))
    message(Sys.time() - t0)
}

# exons ---------------------
myExons <- "../../../../annotationData/gencode.v24.annotation.hg19.middle.exons.pc.light.autosomes.bed"
myExonsTbl <- read_tsv(myExons, col_names = FALSE)
colnames(myExonsTbl) <- c("chr", "start" ,"end", "name", "score" , "strand")

t0 <- Sys.time() # 45 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataForWdith(
        metadata$id[i],
        myExons,
        refgenome = "hg19",
        bin = 50L,
        rm0 = TRUE,
        xmin = 1000L, xmax = 1000L, type = "mf",
        add_heatmap = TRUE,
        verbose = TRUE
    )
    dataForPlot <- cbind(myExonsTbl, gene_type = "protein_coding", dataForPlot) %>% arrange(desc(end - start))

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
        file = paste0("../../../plots/profileHM/exons/byWidth/", metadata$id[i], "_", metadata$short[i], "_exonsByWidth_profile.pdf"),
        width = 8.3, height = 5.8
    )
    plotWidthAndProfile(
        dataForPlot,
        nbin = 9,
        zlims  = zlims,
        raster = TRUE,
        title = gsub("_", " ", metadata$name[i], fixed = TRUE),
        geneTypePalette = colorRampPalette(c("yellow")),
        naColour = "white"
    )
    dev.off()
    message(paste(metadata$id[i], "done!"))
    message(Sys.time() - t0)
}



# tests -----------------


E003data <- addGeneTypeInfo(E003data, refTable)

geneType <- rep("other", nrow(E003data))
geneType[which(E003data$gene_type == "protein_coding")] <- "protein_coding"
geneType[which(E003data$gene_type == "processed_pseudogene")] <- "processed_pseudogene"

E003data2gt <- E003data
E003data2gt$gene_type <- geneType


pdf(file = "../../../plots/E003_heatmap_a5.pdf", width = 8.3, height = 5.8)
plotExpressionAndProfile(
    E003data2gt,
    zlims  = list(c(0, 20)  , c(0, 1)     , c(1, 4)       , c(10, 40)  ),
    raster = TRUE,
    title = "E003 stem cells long string",
    geneTypePalette = colorRampPalette(c("grey", "red", "yellow"))
)
dev.off()






coverageData <- getPlotSetArrayFor(
    "E003",
    myProms,
    refgenome = "hg19",
    bin = 100L,
    rm0 = TRUE,
    xmin = 2500L, xmax = 2500L, type = "pf",
    add_heatmap = TRUE,
    verbose = TRUE
)

roadmapCode <- "E003"

exprTable <-  roadmapExp


exprData <-  select(roadmapExp, gene_id, roadmapCode)





myMat <- combTable[, grep("CpG_sites", colnames(combTable))] %>% as.matrix
layout(rbind(1,2,3), heights = c(1, 0.2, 0.4))

plotTraceAndKeyAndProfile(myMat)


coverageData <- getPlotSetArrayFor(
    metadata$id[i],
    myExons,
    refgenome = "hg19",
    bin = 100L,
    rm0 = TRUE,
    xmin = 2500L, xmax = 2500L, type = "mf",
    add_heatmap = TRUE,
    verbose = TRUE
)

