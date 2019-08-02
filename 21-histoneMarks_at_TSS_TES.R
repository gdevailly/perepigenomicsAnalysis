setwd("/media/gdevailly/SANS TITRE/inra/cascade")

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(parallel)

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/20-functions_for_histoneMarks.R")

# metadata loading, messy ---------------------------------------
metadata <- read_tsv("~/mnt/genotoul_grp/guillaume/cascade/data/wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

myProms <- "~/mnt/genotoul_grp/guillaume/cascade/annotation/gencode.v29.annotation.hg19.middleTSStranscript.light.autosomes.bed"

salmonExp <- read_tsv("~/mnt/genotoul_grp/guillaume/cascade/data/rnaseq/salmon_exp_genes_expressed.tsv")

refTable <- read_tsv(
    "~/mnt/genotoul_grp/guillaume/cascade/annotation/gencode.v29.annotation.hg19.middleTSStranscript.bed",
    col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)
refTable <- filter(refTable, chr %in% paste0("chr", 1:22))

metadata$id[!metadata$id %in% colnames(salmonExp)] # missing E008, E017, E021, E022, and more
metadata <- dplyr::filter(metadata, id %in% colnames(salmonExp))

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
hisFiles <- filter(hisFiles, cellCode %in% colnames(salmonExp))
tss_table <- read_tsv(myProms, col_names = FALSE, progress = FALSE)
colnames(tss_table) <- c("chr", "start", "end", "name", "score", "strand")

adundant_gene_types <- c(
    "protein_coding", "processed_pseudogene", "lincRNA", "antisense",
    "snRNA", "unprocessed_pseudogene", "misc_RNA", "miRNA", "snoRNA",
    "rRNA", "transcribed_unprocessed_pseudogene", "other"
)
refTable$gene_type <- if_else(refTable$gene_type %in% adundant_gene_types, refTable$gene_type, "other")


preffix <- "~/mnt/genotoul/work/projects/cascade/"


# DNAseI -----------------------------
DNAse_md <- inner_join(
    dplyr::filter(hisFiles, ChIP  == "DNase") %>%
        dplyr::select(-name, -ChIP),
    dplyr::filter(hisFiles, ChIP  == "DNaseControl") %>%
        dplyr::select(-name, -ChIP),
    by = "cellCode"
) %>% dplyr::rename(DNAse_file = file.x, Control_file = file.y)

DNAse_md <- filter(DNAse_md, cellCode %in% colnames(salmonExp))

# very unpure function, just for looping actually -------------------
# png versions
plotBetterDNAseDataFor <- function(cellCode, npix_height = 600) {

    dnaseData <- prepareDNAseData(cellCode, tssAnno = tss_table, metadataTable = DNAse_md, expressionTable = salmonExp) %>%
        addGeneTypeInfo(refTable)

    for(i in adundant_gene_types) {
        if(!file.exists(paste0(preffix, "appPlots/tss/dnase1/", i))) {
            dir.create(paste0(preffix, "appPlots/tss/dnase1/", i))
        }
        png(
            file = paste0(
                preffix, "appPlots/tss/dnase1/", i, "/", cellCode, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotMetricAndProfile(
            dplyr::filter(dnaseData, gene_type == i),
            prefix = c("DNAse"  , "Control"),
            labels = c("DNAse", "Control"),
            tints  = c("dodgerblue" , "grey"),
            zlims  = list(c(0, 2) , c(0, 2)),
            title = "",
            withGeneType = FALSE,
            npix_height = npix_height
        )
        dev.off()
    }

    # adding gene type information
    geneType <- rep("other", nrow(dnaseData))
    geneType[which(dnaseData$gene_type == "protein_coding")] <- "protein coding"
    geneType[which(grepl("pseudogene", dnaseData$gene_type, fixed = TRUE))] <- "pseudogene"
    geneType[which(grepl("RNA", dnaseData$gene_type, fixed = TRUE))] <- "RNA gene"
    dnaseData$gene_type <- geneType

    if(!file.exists(paste0(preffix, "appPlots/tss/dnase1/all/"))) {
        dir.create(paste0(preffix, "appPlots/tss/dnase1/all/"))
    }
    png(
        file = paste0(
            preffix, "appPlots/tss/dnase1/all/", cellCode, ".png"
        ),
        width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
    )
    plotMetricAndProfile(
        dnaseData,
        prefix = c("DNAse"  , "Control"),
        labels = c("DNAse", "Control"),
        tints  = c("dodgerblue" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        npix_height = npix_height
    )
    dev.off()

    message(paste(cellCode, "is done!"))
}

# vroom is parralelized
t0 <- Sys.time() # 20 minutes
dummy <- mclapply(
    DNAse_md$cellCode,
    plotBetterDNAseDataFor,
    mc.cores = 2 # big memmory footprint
)
Sys.time() - t0

# danse1 TES ---------------------------
TES <- "~/mnt/genotoul_grp/guillaume/cascade/annotation/gencode.v29.annotation.hg19.middleTES.light.autosomes.bed"
tes_table <- read_tsv(TES, col_names = FALSE, progress = FALSE)
colnames(tes_table) <- c("chr", "start", "end", "name", "score", "strand")
tes_table <- mutate(
    tes_table,
    new_start = if_else(strand == "+", end -1L, start),
    new_end = if_else(strand == "+", end, start + 1L)
) %>%
    dplyr::select(-start, -end) %>%
    dplyr::rename(start = new_start, end = new_end) %>%
    dplyr::select(chr, start, end, name, score, strand)

plotBetterDNAseTESDataFor <- function(cellCode, npix_height = 600) {

    dnaseData <- prepareDNAseData(cellCode, tssAnno = tes_table, metadataTable = DNAse_md, expressionTable = salmonExp) %>%
        addGeneTypeInfo(refTable)

    for(i in adundant_gene_types) {
        if(!file.exists(paste0(preffix, "appPlots/tes/dnase1/", i))) {
            dir.create(paste0(preffix, "appPlots/tes/dnase1/", i))
        }
        png(
            file = paste0(
                preffix, "appPlots/tes/dnase1/", i, "/", cellCode, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13, type = "cairo-png"
        )
        plotMetricAndProfile(
            dplyr::filter(dnaseData, gene_type == i),
            prefix = c("DNAse"  , "Control"),
            labels = c("DNAse", "Control"),
            tints  = c("dodgerblue" , "grey"),
            xlabels = c("-2.5kb", "TES", "+2.5kb"),
            zlims  = list(c(0, 2) , c(0, 2)),
            title = "",
            withGeneType = FALSE,
            npix_height = npix_height
        )
        dev.off()
    }

    # adding gene type information
    geneType <- rep("other", nrow(dnaseData))
    geneType[which(dnaseData$gene_type == "protein_coding")] <- "protein coding"
    geneType[which(grepl("pseudogene", dnaseData$gene_type, fixed = TRUE))] <- "pseudogene"
    geneType[which(grepl("RNA", dnaseData$gene_type, fixed = TRUE))] <- "RNA gene"
    dnaseData$gene_type <- geneType

    if(!file.exists(paste0(preffix, "appPlots/tes/dnase1/all/"))) {
        dir.create(paste0(preffix, "appPlots/tes/dnase1/all/"))
    }
    png(
        file = paste0(
            preffix, "appPlots/tes/dnase1/all/", cellCode, ".png"
        ),
        width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
    )
    plotMetricAndProfile(
        dnaseData,
        prefix = c("DNAse"  , "Control"),
        labels = c("DNAse", "Control"),
        xlabels = c("-2.5kb", "TES", "+2.5kb"),
        tints  = c("dodgerblue" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        npix_height = npix_height
    )
    dev.off()

    message(paste(cellCode, "is done!"))
}

t0 <- Sys.time() # 40 minutes
dummy <- mclapply(
    DNAse_md$cellCode,
    plotBetterDNAseTESDataFor,
    mc.cores = 2 # big memmory footprint
)
Sys.time() - t0

# histone marks -----------------
his_md <- left_join(
    dplyr::filter(hisFiles, !ChIP %in% c("DNase", "DNaseControl", "Input")),
    dplyr::filter(hisFiles, ChIP == "Input") %>%
        dplyr::select(file, cellCode) %>%
        dplyr::rename(input_file = file),
    by = "cellCode"
)

# png for apps ---------------
plotBetterHisModDataForLine <- function(i, npix_height = 600) {

    hisModData <- prepareHisModData(i, tssAnno = tss_table, metadataTable = his_md, expressionTable = salmonExp) %>%
        addGeneTypeInfo(refTable)

    hisModName <- his_md[i, ]$ChIP
    cellID <-  his_md[i, ]$cellCode

    if(!file.exists(paste0(preffix, "appPlots/tss/", hisModName))) {
        dir.create(paste0(preffix, "appPlots/tss/", hisModName))
    }

    for(j in adundant_gene_types) {
        if(!file.exists(paste0(preffix, "appPlots/tss/", hisModName, "/", j))) {
            dir.create(paste0(preffix, "appPlots/tss/", hisModName, "/", j))
        }
        png(
            file = paste0(
                preffix, "appPlots/tss/", hisModName, "/", j, "/", cellID, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotMetricAndProfile(
            dplyr::filter(hisModData, gene_type == j),
            prefix = c("HisMod"  , "Input"),
            labels = c(hisModName, "Input"),
            tints  = c("sienna1" , "grey"),
            zlims  = list(c(0, 2) , c(0, 2)),
            title = "",
            withGeneType = FALSE,
            npix_height = npix_height
        )
        dev.off()
    }

    # adding gene type information
    geneType <- rep("other", nrow(hisModData))
    geneType[which(hisModData$gene_type == "protein_coding")] <- "protein coding"
    geneType[which(grepl("pseudogene", hisModData$gene_type, fixed = TRUE))] <- "pseudogene"
    geneType[which(grepl("RNA", hisModData$gene_type, fixed = TRUE))] <- "RNA gene"
    hisModData$gene_type <- geneType

    if(!file.exists(paste0(preffix, "appPlots/tss/", hisModName, "/all/"))) {
        dir.create(paste0(preffix, "appPlots/tss/", hisModName, "/all/"))
    }

    png(
        file = paste0(
            preffix, "appPlots/tss/", hisModName, "/all/", cellID, ".png"
        ),
        width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
    )
    plotMetricAndProfile(
        hisModData,
        prefix = c("HisMod"  , "Input"),
        labels = c(hisModName, "Input"),
        tints  = c("sienna1" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        npix_height = npix_height
    )
    dev.off()

    message(paste(his_md[i, ]$name, "is done!"))
}

t0 <- Sys.time() # 9h... c'est long
dummy <- lapply(
    137:242,
    plotBetterHisModDataForLine
)
Sys.time() - t0

t0 <- Sys.time() # 9h... c'est long
dummy <- lapply(
    seq_len(nrow(his_md)),
    plotBetterHisModDataForLine
)
Sys.time() - t0


# tes -----------------------------
# histones tes ---------------------
TES <- "~/mnt/genotoul_grp/guillaume/cascade/annotation/gencode.v29.annotation.hg19.middleTES.light.autosomes.bed"
tes_table <- read_tsv(TES, col_names = FALSE, progress = FALSE)
colnames(tes_table) <- c("chr", "start", "end", "name", "score", "strand")
tes_table <- mutate(
    tes_table,
    new_start = if_else(strand == "+", end -1L, start),
    new_end = if_else(strand == "+", end, start + 1L)
) %>%
    dplyr::select(-start, -end) %>%
    dplyr::rename(start = new_start, end = new_end) %>%
    dplyr::select(chr, start, end, name, score, strand)

# png for apps ---------------
plotBetterHisModDataTesForLine <- function(i, npix_height = 600) {

    hisModData <- prepareHisModData(i, tssAnno = tes_table, metadataTable = his_md, expressionTable = salmonExp) %>%
        addGeneTypeInfo(refTable)

    hisModName <- his_md[i, ]$ChIP
    cellID <-  his_md[i, ]$cellCode

    if(!file.exists(paste0(preffix, "appPlots/tes/", hisModName))) {
        dir.create(paste0(preffix, "appPlots/tes/", hisModName))
    }

    for(j in adundant_gene_types) {
        if(!file.exists(paste0(preffix, "appPlots/tes/", hisModName, "/", j))) {
            dir.create(paste0(preffix, "appPlots/tes/", hisModName, "/", j))
        }
        png(
            file = paste0(
                preffix, "appPlots/tes/", hisModName, "/", j, "/", cellID, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotMetricAndProfile(
            dplyr::filter(hisModData, gene_type == j),
            prefix = c("HisMod"  , "Input"),
            labels = c(hisModName, "Input"),
            tints  = c("sienna1" , "grey"),
            zlims  = list(c(0, 2) , c(0, 2)),
            title = "",
            withGeneType = FALSE,
            npix_height = npix_height,
            xlabels = c("-2.5kb", "TES", "+2.5kb")
        )
        dev.off()
    }

    # adding gene type information
    geneType <- rep("other", nrow(hisModData))
    geneType[which(hisModData$gene_type == "protein_coding")] <- "protein coding"
    geneType[which(grepl("pseudogene", hisModData$gene_type, fixed = TRUE))] <- "pseudogene"
    geneType[which(grepl("RNA", hisModData$gene_type, fixed = TRUE))] <- "RNA gene"
    hisModData$gene_type <- geneType

    if(!file.exists(paste0(preffix, "appPlots/tes/", hisModName, "/all/"))) {
        dir.create(paste0(preffix, "appPlots/tes/", hisModName, "/all/"))
    }

    png(
        file = paste0(
            preffix, "appPlots/tes/", hisModName, "/all/", cellID, ".png"
        ),
        width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
    )
    plotMetricAndProfile(
        hisModData,
        prefix = c("HisMod"  , "Input"),
        labels = c(hisModName, "Input"),
        tints  = c("sienna1" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        npix_height = npix_height,
        xlabels = c("-2.5kb", "TES", "+2.5kb")
    )
    dev.off()

    rm(hisModData)
    gc()

    message(paste(his_md[i, ]$name, "is done!"))
}

t0 <- Sys.time() # 6h... c'est long
dummy <- lapply(
    seq_len(nrow(his_md)),
    plotBetterHisModDataTesForLine
)
Sys.time() - t0



t0 <- Sys.time() # 6h... c'est long
dummy <- lapply(
    136:242,
    plotBetterHisModDataTesForLine
)
Sys.time() - t0

