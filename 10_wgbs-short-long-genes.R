setwd("~/mnt/genotoul_grp/guillaume/cascade")


library(dplyr)
library(plotrix)
library(seqplots)
library(readr)
library(svglite)
library(parallel)

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")

metadata <- read_tsv("data/wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

myProms <- "annotation/gencode.v29.annotation.hg19.middleTSStranscript.light.autosomes.bed"
salmonExp <- read_tsv("data/rnaseq/salmon_exp_genes_expressed.tsv")

refTable <- read_tsv(
    "annotation/gencode.v29.annotation.hg19.middleTSStranscript.bed",
    col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)
refTable <- dplyr::filter(refTable, chr %in% paste0("chr", 1:22))
metadata$id[!metadata$id %in% colnames(salmonExp)] # missing E008, E017, E021, E022, and others
metadata <- dplyr::filter(metadata, id %in% colnames(salmonExp))
adundant_gene_types <- c(
    "protein_coding", "processed_pseudogene", "lincRNA", "antisense",
    "snRNA", "unprocessed_pseudogene", "misc_RNA", "miRNA", "snoRNA",
    "rRNA", "transcribed_unprocessed_pseudogene", "other"
)
refTable$gene_type <- if_else(refTable$gene_type %in% adundant_gene_types, refTable$gene_type, "other")
preffix <- "~/mnt/genotoul/work/projects/cascade/"

# long and short genes ------------
summary(refTable$end - refTable$start)
refTable <- mutate(
    refTable,
    gene_length = end - start
)
refTable <- mutate(refTable, length_type = case_when(
    refTable$gene_length <= 1000 ~ "short",
    refTable$gene_length <= 3000 ~ "intermediate",
    refTable$gene_length > 3000 ~ "long",
    TRUE ~ NA_character_
))
table(refTable$length_type, useNA = "ifany")
#intermediate         long        short
# 9335        25954        20112


plotLengthTypeWgbsDataFor <- function(i, npix_height = 600) {

    gc()

    dataForPlot <- extractAndPrepareDataFor(
        metadata$id[i],
        myProms,
        salmonExp,
        refgenome = "hg19",
        bin = 100L,
        rm0 = TRUE,
        xmin = 2500L, xmax = 2500L, type = "pf",
        add_heatmap = TRUE,
        verbose = FALSE
    ) %>% addLengthTypeInfo(refTable)

    # adding gene type information
    geneType <- rep("other", nrow(dataForPlot))
    geneType[which(dataForPlot$gene_type == "protein_coding")] <- "protein coding"
    geneType[which(grepl("pseudogene", dataForPlot$gene_type, fixed = TRUE))] <- "pseudogene"
    geneType[which(grepl("RNA", dataForPlot$gene_type, fixed = TRUE))] <- "RNA gene"
    dataForPlot$gene_type <- geneType

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

    for(j in c("short", "intermediate", "long")) {
        if(!file.exists(paste0(preffix, "appPlots/tss/wgbs/", j))) {
            dir.create(paste0(preffix, "appPlots/tss/wgbs/", j))
        }
        png(
            file = paste0(
                preffix, "appPlots/tss/wgbs/", j, "/", metadata$id[i], ".png"
            ),
            width = 8.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotMetricAndProfile(
            dplyr::filter(dataForPlot, length_type == j),
            zlims  = zlims,
            raster = TRUE,
            title = "",
            tints  = c("red", "blue", "purple", "grey"),
            withGeneType = TRUE,
            npix_height = npix_height
        )
        dev.off()
    }

    message(paste(metadata$id[i], "is done!"))
}

t0 <- Sys.time() # 13 minutes
dummy <- mclapply(
    seq_len(nrow(metadata)),
    plotLengthTypeWgbsDataFor,
    mc.cores = 2
)
Sys.time() - t0


# TES ----------------------
TES <- "annotation/gencode.v29.annotation.hg19.middleTES.light.autosomes.bed"

plotBetterWgbsDataTESFor <- function(i, npix_height = 600) {

    dataForPlot <- extractAndPrepareDataFor(
        metadata$id[i],
        TES,
        salmonExp,
        refgenome = "hg19",
        bin = 100L,
        rm0 = TRUE,
        xmin = 2500L, xmax = 2500L, type = "ef",
        add_heatmap = TRUE,
        verbose = FALSE
    ) %>% addLengthTypeInfo(refTable)

    # adding gene type information
    geneType <- rep("other", nrow(dataForPlot))
    geneType[which(dataForPlot$gene_type == "protein_coding")] <- "protein coding"
    geneType[which(grepl("pseudogene", dataForPlot$gene_type, fixed = TRUE))] <- "pseudogene"
    geneType[which(grepl("RNA", dataForPlot$gene_type, fixed = TRUE))] <- "RNA gene"
    dataForPlot$gene_type <- geneType

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

    for(j in c("short", "intermediate", "long")) {
        if(!file.exists(paste0(preffix, "appPlots/tes/wgbs/", j))) {
            dir.create(paste0(preffix, "appPlots/tes/wgbs/", j))
        }
        png(
            file = paste0(
                preffix, "appPlots/tes/wgbs/", j, "/", metadata$id[i], ".png"
            ),
            width = 8.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotMetricAndProfile(
            dplyr::filter(dataForPlot, length_type == j),
            zlims  = zlims,
            raster = TRUE,
            title = "",
            tints  = c("red", "blue", "purple", "grey"),
            xlabels = c("-2.5kb", "TTS", "+2.5kb"),
            withGeneType = TRUE,
            npix_height = npix_height
        )
        dev.off()
    }

    message(paste(metadata$id[i], "is done!"))
}

t0 <- Sys.time() # 12 minutes
dummy <- mclapply(
    seq_len(nrow(metadata)),
    plotBetterWgbsDataTESFor,
    mc.cores = 2
)
Sys.time() - t0

