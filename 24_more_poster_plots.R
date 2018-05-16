# DNAme at TES ---------------
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
TES <- "../../../../annotationData/gencode.v24.annotation.hg19.middleTES.light.autosomes.bed"

i <- 1 #hESC
dataForPlot <- extractAndPrepareDataFor(
    metadata$id[i],
    TES,
    roadmapExp,
    refgenome = "hg19",
    bin = 100L,
    rm0 = TRUE,
    xmin = 2000L, xmax = 2000L, type = "ef",
    add_heatmap = TRUE,
    verbose = TRUE
)
dataForPlot <- addGeneTypeInfo(dataForPlot, refTable)

plotWidth <- 1.25
png(
    file = paste0("../../../plots/poster/hESC_tts_meDNA.png"),
    width = 4.1*plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
plotProfileOnly(
    dataForPlot,
    zlims  = zlims <- list(c(0, 15), c(0,1), c(1, 4), c(10, 50)),
    raster = FALSE,
    title = "DNA methylation",
    naColour = "grey",
    tints  = c("red", "blue", "purple", "grey"),
    xlabels = c("-2kb", "TTS", "+2kb"),
    npix_height = 400
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

plotWidth <- 1.25
png(
    file = paste0("../../../plots/poster/hESC_tts_meDNA_protCod.png"),
    width = 4.1*plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
plotProfileOnly(
    dplyr::filter(dataForPlot, gene_type == "protein_coding"),
    zlims  = zlims <- list(c(0, 15), c(0,1), c(1, 4), c(10, 50)),
    raster = FALSE,
    title = "DNA methylation",
    naColour = "grey",
    tints  = c("red", "blue", "purple", "grey"),
    xlabels = c("-2kb", "TTS", "+2kb"),
    npix_height = 400
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# h3K36me3 (and H3K4me3 as a negative control) at TES ---------------------

setwd("/groups2/joshi_grp/guillaume/cascade/data/")

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(parallel)

source("../Rscripts/6-plotingFunctions.R")
source("../Rscripts/20-functions_for_histoneMarks.R")

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

TES <- "../../annotationData/gencode.v24.annotation.hg19.middleTES.light.autosomes.bed"
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
his_md <- left_join(
    dplyr::filter(hisFiles, !ChIP %in% c("DNase", "DNaseControl", "Input")),
    dplyr::filter(hisFiles, ChIP == "Input") %>%
        dplyr::select(file, cellCode) %>%
        dplyr::rename(input_file = file),
    by = "cellCode"
)

hisH3k36me3 <- filter(
    his_md, cellCode == "E003", grepl("H3K36me3", ChIP)
)
hisH3k4me3 <- filter(
    his_md, cellCode == "E003", grepl("H3K4me3", ChIP)
)

t0 <- Sys.time() # 6 minutes
signalMat <- importReadAlignAsGranges(
    hisH3k36me3$file
) %>% extractFeatureScoreMatrix(annoTable = tes_table, up = 2000, down = 2000)
signalDf <- data.frame(
    gene_id = rownames(signalMat),
    signalMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(signalDf) <- c("gene_id", paste0(hisH3k36me3$ChIP, colnames(signalDf)[-1]))
signalDf <- inner_join(
    select_(roadmapExp, "gene_id",  hisH3k36me3$cellCode),
    signalDf,
    by = "gene_id"
) %>%
    rename_(exp = hisH3k36me3$cellCode) %>%
    arrange(desc(exp))

signalMat2 <- importReadAlignAsGranges(
    hisH3k4me3$file
) %>% extractFeatureScoreMatrix(annoTable = tes_table, up = 2000, down = 2000)
signalDf2 <- data.frame(
    gene_id = rownames(signalMat2),
    signalMat2
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(signalDf2) <- c("gene_id", paste0(hisH3k4me3$ChIP, colnames(signalDf2)[-1]))
signalDf2 <- inner_join(
    select_(roadmapExp, "gene_id",  hisH3k4me3$cellCode),
    signalDf2,
    by = "gene_id"
) %>%
    rename_(exp = hisH3k4me3$cellCode) %>%
    arrange(desc(exp))

inputMat <- importReadAlignAsGranges(
    hisH3k36me3$input_file
) %>% extractFeatureScoreMatrix(annoTable = tes_table, up = 2000, down = 2000)
inputDf <- data.frame(
    gene_id = rownames(inputMat),
    inputMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(inputDf) <- c("gene_id", paste0("Input", colnames(inputDf)[-1]))
inputDf <- inner_join(
    select_(roadmapExp, "gene_id",  hisH3k36me3$cellCode),
    inputDf,
    by = "gene_id"
) %>%
    rename_(exp = hisH3k36me3$cellCode) %>%
    arrange(desc(exp))

myhisH3k36me3TabClean <- cbind(
    signalDf,
    dplyr::select(signalDf2, -gene_id, - exp),
    dplyr::select(inputDf  , -gene_id, - exp)
)
Sys.time() - t0

myprefix <- c("H3K4me3", "H3K36me3", "Input")
plotWidth <- 1.25
png(
    file = "../plots/poster/hESC_tts_H3K36me3.png",
    width = (length(myprefix)+0.1)*plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
plotProfileOnly(
    myhisH3k36me3TabClean,
    zlims  = list(c(0, 1), c(0, 1), c(0, 1)),
    raster = TRUE,
    title = "Histone methylation",
    prefix = myprefix,
    labels = c("H3K4me3", "H3K36me3", "ChIP input"),
    tints  = c("sienna1", "sienna1", "grey"),
    naColour = "grey",
    xlabels = c("-2kb", "TTS", "+2kb"),
    npix_height = 500
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# DNAme at exons ---------------
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

i <- 1
dataForPlot <- extractAndPrepareDataForExons(
    metadata$id[i],
    exonInfoPath,
    innerExonRatio,
    refgenome = "hg19",
    bin = 50L,
    rm0 = TRUE,
    xmin = 500L, xmax = 500L, type = "pf",
    add_heatmap = TRUE,
    verbose = TRUE
)

zlims <- list(c(0, 15), c(0, 1), c(1, 5), c(10, 50))
bins <- numeric(21125)
bins[which(dataForPlot$exp >= 2)] <- 1
bins[which(dataForPlot$exp < 2 & dataForPlot$exp >= 1)] <- 2
bins[which(dataForPlot$exp < 1 & dataForPlot$exp >= 1/2)] <- 3
bins[which(dataForPlot$exp < 1/2)] <- 4
plotWidth <- 1.25

png(
    file = "../../../plots/poster/hESC_exonsRatio_DNAme.png",
    width = (4)*plotWidth+0.1+1.7, height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
plotMetricAndProfile(
    dataForPlot,
    bins = bins,
    zlims  = zlims,
    raster = TRUE,
    title = "DNA methylation",
    geneTypePalette = colorRampPalette(c("yellow")),
    withGeneType = FALSE,
    expTransf = function(x) log2(x+0.01),
    xlabels = c("-500", "Exon st.", "+500"),
    expLim = c(-3, 3),
    exp.text = "log2(exon/gene)",
    reversedZOrder = TRUE,
    tints  = c("red", "blue", "purple", "grey"),
    npix_height = 400,
    abline.v = 0
)
dev.off()

# h3k36me3 at exons fpkm --------------------------
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

# build histone marks table
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

# DNAseI
DNAse_md <- inner_join(
    dplyr::filter(hisFiles, ChIP  == "DNase") %>%
        dplyr::select(-name, -ChIP),
    dplyr::filter(hisFiles, ChIP  == "DNaseControl") %>%
        dplyr::select(-name, -ChIP),
    by = "cellCode"
) %>% dplyr::rename(DNAse_file = file.x, Control_file = file.y)


exonInfoPath <- "inner_exon.bed"
exonTable <- read_tsv(exonInfoPath, col_names = FALSE, progress = FALSE)
colnames(exonTable) <- c("chr", "start", "end", "name", "score", "strand")
innerExonQuant <- dplyr::select(innerExonQuant, -name) %>%
    dplyr::rename(name = exon_location)

hisH3k36me3 <- filter(
    his_md, cellCode == "E003", grepl("H3K36me3", ChIP)
)
hisH3k4me3 <- filter(
    his_md, cellCode == "E003", grepl("H3K4me3", ChIP)
)

t0 <- Sys.time() # 6 minutes
signalMat <- importReadAlignAsGranges(
    hisH3k36me3$file
) %>% extractFeatureScoreMatrix(annoTable = exonTable, up = 500, down = 500, freq = 50)
signalDf <- data.frame(
    gene_id = rownames(signalMat),
    signalMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(signalDf) <- c("gene_id", paste0(hisH3k36me3$ChIP, colnames(signalDf)[-1]))
signalDf <- inner_join(
    select_(innerExonQuant, "name",  hisH3k36me3$cellCode),
    signalDf,
    by = c("name" = "gene_id")
) %>%
    rename_(exp = hisH3k36me3$cellCode) %>%
    arrange(desc(exp))

signalMat2 <- importReadAlignAsGranges(
    hisH3k4me3$file
) %>% extractFeatureScoreMatrix(annoTable = exonTable, up = 2000, down = 2000, freq = 50)
signalDf2 <- data.frame(
    gene_id = rownames(signalMat2),
    signalMat2
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(signalDf2) <- c("gene_id", paste0(hisH3k4me3$ChIP, colnames(signalDf2)[-1]))
signalDf2 <- inner_join(
    select_(innerExonQuant, "name",  hisH3k4me3$cellCode),
    signalDf2,
    by = c("name" = "gene_id")
) %>%
    rename_(exp = hisH3k4me3$cellCode) %>%
    arrange(desc(exp))

inputMat <- importReadAlignAsGranges(
    hisH3k36me3$input_file
) %>% extractFeatureScoreMatrix(annoTable = exonTable, up = 2000, down = 2000, freq = 50)
inputDf <- data.frame(
    gene_id = rownames(inputMat),
    inputMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(inputDf) <- c("gene_id", paste0("Input", colnames(inputDf)[-1]))
inputDf <- inner_join(
    select_(innerExonQuant, "name",  hisH3k36me3$cellCode),
    inputDf,
    by = c("name" = "gene_id")
) %>%
    rename_(exp = hisH3k36me3$cellCode) %>%
    arrange(desc(exp))

myhisH3k36me3TabClean <- cbind(
    signalDf,
    dplyr::select(signalDf2, -name, - exp),
    dplyr::select(inputDf  , -name, - exp)
)
Sys.time() - t0

plotWidth <- 1.25
png(
    file = "../plots/poster/hESC_exons_fpkm_H3K36me3.png",
    width = (3)*plotWidth+0.1+1.7, height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
plotMetricAndProfile(
    myhisH3k36me3TabClean,
    zlims  = list(c(0,2), c(0,2), c(0,1)),
    raster = TRUE,
    title = "Histone methylation",
    withGeneType = FALSE,
    xlabels = c("-500", "Exon st.", "+500"),
    expLim = c(0, 3.5),
    exp.text = "log10(exon FPKM+1)",
    reversedZOrder = FALSE,
    prefix = c("H3K4me3", "H3K36me3", "Input"),
    labels = c("H3K4me3", "H3K36me3", "ChIP input"),
    tints  = c("sienna1", "sienna1", "grey"),
    npix_height = 500
)
dev.off()


