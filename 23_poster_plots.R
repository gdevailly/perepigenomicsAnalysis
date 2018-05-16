# DNAme at TSS, prepare data -----------------
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

i <- 1 #hESC
dataForPlot <- extractAndPrepareDataFor(
    metadata$id[i],
    myProms,
    roadmapExp,
    refgenome = "hg19",
    bin = 100L,
    rm0 = TRUE,
    xmin = 2000L, xmax = 2000L, type = "pf",
    add_heatmap = TRUE,
    verbose = TRUE
)
dataForPlot <- addGeneTypeInfo(dataForPlot, refTable)

geneType <- rep("other", nrow(dataForPlot))
geneType[which(dataForPlot$gene_type == "protein_coding")] <- "protein coding"
geneType[which(grepl("pseudogene", dataForPlot$gene_type, fixed = TRUE))] <- "pseudogene"
geneType[which(grepl("RNA", dataForPlot$gene_type, fixed = TRUE))] <- "RNA gene"
dataForPlot$gene_type <- geneType
myDF <- dataForPlot

# expression plots ---------------------
png(
    file = paste0("../../../plots/poster/meDNA_", metadata$id[i], "_", metadata$short[i], "_tss_geneExp_BIG.png"),
    width = 1.7, height = 6, units = "in", res = 1200, pointsize = 13, bg = "transparent"#, type = "Xlib"#, antialias = "none"
)
layout(
    matrix(c(1, 4, 2, 3, 4, 2, 5, 4, 6), nrow = 3),
    widths = c(1, 0.2, 0.2),
    c(0.65, 0.15, 0.30)
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
# expression
plotExpression2(myDF, expLim = c(0, 3.5))
# gene type
par(mar = c(0, 0.8, 4, 0))
geneTypePalette <- colorRampPalette(c("grey", "blue", "yellow", "red"))
geneType <- as.numeric(factor(myDF$gene_type))
image(t(rev(geneType)), col = geneTypePalette(length(unique(myDF$gene_type))), axes = FALSE,
      breaks = seq(0.5, length(unique(myDF$gene_type)) + 0.5, by = 1))
box()
mtext(side = 3, "gene\ntype", line = 0.5,  cex = 0.8, las = 2)
par(mar = c(0, 0, 0, 0))
plot.new()
ml <- legend("topright", legend = rep("", length(unique(myDF$gene_type))), fill = geneTypePalette(length(unique(myDF$gene_type))), bty = "n",
             title = "gene type:", title.adj = -5)
text(ml$text$x - 0.08, ml$text$y - 0.005, levels(factor(myDF$gene_type)), pos = 2)
# bining
par(mar = c(0, 0.8, 4, 0))
image(t(rev(ntile(seq_len(nrow(myDF)), 5))), col = colorRampPalette(c("magenta", "black", "green"))(5), axes = FALSE)
axis(4, at = seq(0, 1, length.out = 5 + 1), labels = NA)
box()
mtext(side = 3, "bin", line = 0.5,  cex = 0.8, las = 2)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# DNAme plots ------------------
plotWidth <- 1.25
png(
    file = paste0("../../../plots/poster/", metadata$id[i], "_", metadata$short[i], "_tss_meDNA.png"),
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
    xlabels = c("-2kb", "TSS", "+2kb"),
    npix_height = 400
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# new R session ---------------------------------
# DNAseI ----------------------------------
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

tss_table <- read_tsv(myProms, col_names = FALSE, progress = FALSE)
colnames(tss_table) <- c("chr", "start", "end", "name", "score", "strand")

# DNAseI -----------------------------
DNAse_md <- inner_join(
    dplyr::filter(hisFiles, ChIP  == "DNase") %>%
        dplyr::select(-name, -ChIP),
    dplyr::filter(hisFiles, ChIP  == "DNaseControl") %>%
        dplyr::select(-name, -ChIP),
    by = "cellCode"
) %>% dplyr::rename(DNAse_file = file.x, Control_file = file.y)

cellCode <- "E003"
dnaseData <- prepareDNAseData(cellCode, tssAnno = tss_table, metadataTable = DNAse_md, expressionTable = roadmapExp, up = 2000, down = 2000) %>%
    addGeneTypeInfo(refTable)

# adding gene type information
geneType <- rep("other", nrow(dnaseData))
geneType[which(dnaseData$gene_type == "protein_coding")] <- "protein coding"
geneType[which(grepl("pseudogene", dnaseData$gene_type, fixed = TRUE))] <- "pseudogene"
geneType[which(grepl("RNA", dnaseData$gene_type, fixed = TRUE))] <- "RNA gene"
dnaseData$gene_type <- geneType

plotWidth <- 1.25
png(
    file = "../plots/poster/hESC_tss_DNAse1.png",
    width = 2.1*plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
plotProfileOnly(
    dnaseData,
    zlims  = list(c(0, 3) , c(0, 1)),
    raster = TRUE,
    title = "DNAseI",
    prefix = c("DNAse"  , "Control"),
    labels = c("DNAseI", "Control"),
    tints  = c("darkgreen" , "grey"),
    naColour = "grey",
    xlabels = c("-2kb", "TSS", "+2kb"),
    npix_height = 500
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# histome marks -------
# acetylations --------
hisAc <- filter(
    his_md, cellCode == "E003", grepl("ac$", ChIP)
)

t0 <- Sys.time() # 3 minutes
signalMats <- mclapply(
    seq_len(nrow(hisAc)),
    function(i) {
        signalMat <- importReadAlignAsGranges(
            hisAc[i,]$file
        ) %>% extractFeatureScoreMatrix(annoTable = tss_table, up = 2000, down = 2000)
        signalDf <- data.frame(
            gene_id = rownames(signalMat),
            signalMat
        ) %>% as_data_frame %>%
            mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
        colnames(signalDf) <- c("gene_id", paste0(hisAc[i,]$ChIP, colnames(signalDf)[-1]))
        return(
            inner_join(
                select_(roadmapExp, "gene_id",  hisAc[i,]$cellCode),
                signalDf,
                by = "gene_id"
            ) %>%
                rename_(exp = hisAc[i,]$cellCode) %>%
                arrange(desc(exp))
        )
    },
    mc.cores = 20
)
Sys.time() - t0

myHisAcTab <- signalMats[[1]]
signalMats <- signalMats[-1]
myHisAcTabClean <- cbind(
    myHisAcTab,
    do.call(cbind, lapply(signalMats, function(x) dplyr::select(x, -gene_id, - exp)))
)

myprefix <- hisAc$ChIP
myprefix <- c(
    "H2AK5ac", "H2BK5ac", "H2BK12ac", "H2BK15ac", "H2BK20ac", "H2BK120ac",
    "H3K4ac", "H3K9ac", "H3K14ac", "H3K18ac", "H3K23ac", "H3K27ac", "H3K56ac",
    "H4K5ac", "H4K8ac", "H4K91ac"
)
plotWidth <- 1.25
zlims <- map(myprefix, function(x) return(c(0, 1)))

png(
    file = "../plots/poster/hESC_tss_HisAc.png",
    width = (length(myprefix)+0.1)*plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
plotProfileOnly(
    myHisAcTabClean,
    zlims  = zlims,
    raster = TRUE,
    title = "Histone Acetylation",
    prefix = myprefix,
    labels = myprefix,
    tints  = rep("darkorange", length(myprefix)),
    naColour = "grey",
    xlabels = c("-2kb", "TSS", "+2kb"),
    npix_height = 500
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# metylations ----------------------
hisMe <- filter(
    his_md, cellCode == "E003", grepl("me", ChIP)
)

t0 <- Sys.time() # 6 minutes
signalMats <- mclapply(
    seq_len(nrow(hisMe)),
    function(i) {
        signalMat <- importReadAlignAsGranges(
            hisMe[i,]$file
        ) %>% extractFeatureScoreMatrix(annoTable = tss_table, up = 2000, down = 2000)
        signalDf <- data.frame(
            gene_id = rownames(signalMat),
            signalMat
        ) %>% as_data_frame %>%
            mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
        colnames(signalDf) <- c("gene_id", paste0(hisMe[i,]$ChIP, colnames(signalDf)[-1]))
        return(
            inner_join(
                select_(roadmapExp, "gene_id",  hisMe[i,]$cellCode),
                signalDf,
                by = "gene_id"
            ) %>%
                rename_(exp = hisMe[i,]$cellCode) %>%
                arrange(desc(exp))
        )
    },
    mc.cores = 20
)
Sys.time() - t0

myHisMeTab <- signalMats[[1]]
signalMats <- signalMats[-1]
myHisMeTabClean <- cbind(
    myHisMeTab,
    do.call(cbind, lapply(signalMats, function(x) dplyr::select(x, -gene_id, - exp)))
)

myprefix <- hisMe$ChIP
myprefix <- c(
    "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3", "H3K23me2", "H3K27me3", "H3K36me3",
    "H3K79me1", "H3K79me2", "H4K20me1"
)
plotWidth <- 1.25
zlims <- map(myprefix, function(x) return(c(0, 1)))
zlims[[2]] <- c(0, 4)
zlims[[3]] <- c(0, 4)
png(
    file = "../plots/poster/hESC_tss_HisMe.png",
    width = (length(myprefix)+0.1)*plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
plotProfileOnly(
    myHisMeTabClean,
    zlims  = zlims,
    raster = TRUE,
    title = "Histone Methylation",
    prefix = myprefix,
    labels = myprefix,
    tints  = rep("sienna1", length(myprefix)),
    naColour = "grey",
    xlabels = c("-2kb", "TSS", "+2kb"),
    npix_height = 500
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# H2AZ + Input -----------------
hisH2AZ <- filter(
    his_md, cellCode == "E003", grepl("H2A.Z", ChIP)
)

t0 <- Sys.time() # 6 minutes
signalMat <- importReadAlignAsGranges(
    hisH2AZ$file
) %>% extractFeatureScoreMatrix(annoTable = tss_table, up = 2000, down = 2000)
signalDf <- data.frame(
    gene_id = rownames(signalMat),
    signalMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(signalDf) <- c("gene_id", paste0(hisH2AZ$ChIP, colnames(signalDf)[-1]))
signalDf <- inner_join(
    select_(roadmapExp, "gene_id",  hisH2AZ$cellCode),
    signalDf,
    by = "gene_id"
) %>%
    rename_(exp = hisH2AZ$cellCode) %>%
    arrange(desc(exp))

inputMat <- importReadAlignAsGranges(
    hisH2AZ$input_file
) %>% extractFeatureScoreMatrix(annoTable = tss_table, up = 2000, down = 2000)
inputDf <- data.frame(
    gene_id = rownames(inputMat),
    inputMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(inputDf) <- c("gene_id", paste0("Input", colnames(inputDf)[-1]))
inputDf <- inner_join(
    select_(roadmapExp, "gene_id",  hisH2AZ$cellCode),
    inputDf,
    by = "gene_id"
) %>%
    rename_(exp = hisH2AZ$cellCode) %>%
    arrange(desc(exp))
myHisH2AZTabClean <- cbind(
    signalDf,
    dplyr::select(inputDf, -gene_id, - exp)
)
Sys.time() - t0

myprefix <- c("H2A.Z", "Input")
plotWidth <- 1.25
png(
    file = "../plots/poster/hESC_tss_H2AZ.png",
    width = (length(myprefix)+0.1)*plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
)
oldMar <- par()$mar
oldMgp <- par()$mgp
par(mgp = c(1.4, 0.5, 0))
plotProfileOnly(
    myHisH2AZTabClean,
    zlims  = list(c(0, 2), c(0, 1)),
    raster = TRUE,
    title = "",
    prefix = myprefix,
    labels = c("H2A.Z", "ChIP input"),
    tints  = c("dodgerblue", "grey"),
    naColour = "grey",
    xlabels = c("-2kb", "TSS", "+2kb"),
    npix_height = 500
)
par(mar = oldMar, mgp = oldMgp)
dev.off()

# h3K27me3 by gene type --------------------
hish3k27me3 <- filter(
    his_md, cellCode == "E003", grepl("H3K27me3", ChIP)
)
signalMat <- importReadAlignAsGranges(
    hish3k27me3$file
) %>% extractFeatureScoreMatrix(annoTable = tss_table, up = 2000, down = 2000)
signalDf <- data.frame(
    gene_id = rownames(signalMat),
    signalMat
) %>% as_data_frame %>%
    mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
colnames(signalDf) <- c("gene_id", paste0(hish3k27me3$ChIP, colnames(signalDf)[-1]))
signalDf <- inner_join(
    select_(roadmapExp, "gene_id",  hish3k27me3$cellCode),
    signalDf,
    by = "gene_id"
) %>%
    rename_(exp = hish3k27me3$cellCode) %>%
    arrange(desc(exp))
signalDf <- addGeneTypeInfo(signalDf, refTable)
table(signalDf$gene_type) %>% sort %>% rev

myTables <- map(
    c("protein_coding", "processed_pseudogene", "lincRNA", "antisense", "snRNA", "unprocessed_pseudogene",
      "misc_RNA", "miRNA", "snoRNA", "rRNA"),
    ~dplyr::filter(signalDf, gene_type == .x)
)
names(myTables) <- c("protein_coding", "processed_pseudogene", "lincRNA", "antisense", "snRNA", "unprocessed_pseudogene",
                     "misc_RNA", "miRNA", "snoRNA", "rRNA")

lapply(
    names(myTables),
    function(x) {
        png(
            file = paste0("../plots/poster/hESC_tss_H3K27_", x, ".png"),
            width = 1.7+plotWidth , height = 6, pointsize = 13, bg = "transparent", units = "in", res = 600
        )
        plotMetricAndProfile(
            myTables[[x]],
            nbin = 5,
            bins = NULL, # nbin is ignored if bin is specified. Bins should be a numeric vector of length nrow(myDF)
            prefix = c("H3K27me3"),
            labels = c("H3K27me3"),
            tints  = c("sienna1"        , "blue"      , "purple"      , "green"   ),
            zlims  = list(c(0, 1)),
            naColour = "grey",
            withGeneType = FALSE,
            title = x,
            expTransf = function(x) log10(x+1),
            expLim = c(0, 3.5),
            exp.text = "log10(FPKM+1)",
            xlabels = c("-2kb", "TSS", "+2kb"),
            npix_height = 500
        )
        par(mar = oldMar, mgp = oldMgp)
        dev.off()
    }
) %>% invisible

