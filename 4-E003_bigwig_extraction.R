setwd("/groups2/joshi_grp/guillaume/cascade/data/wgbs/roadmap")
library(seqplots)

myBW <- c(
    "E003_WGBS_FractionalMethylation.bigwig",
    "E003_CpG_sites_hg19.bw",
    "E003_mCpG_density_hg19.bw",
    "E003_mCpG_ratio_hg19.bw",
    "E003_WGBS_ReadCoverage.bigwig"
)
myProms <- "../../../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.light.autosomes.bed"

t0 <- Sys.time()
E003_psa_no0 <- getPlotSetArray(
    myBW,
    myProms,
    refgenome = "hg19",
    bin = 100L,
    rm0 = TRUE,
    xmin = 2500L, xmax = 2500L, type = "pf",
    add_heatmap = TRUE,
    verbose = TRUE
)
Sys.time() - t0 # 1 min

t0 <- Sys.time()
E003_psa_w0 <- getPlotSetArray(
    myBW,
    myProms,
    refgenome = "hg19",
    bin = 100L,
    rm0 = FALSE,
    xmin = 2500L, xmax = 2500L, type = "pf",
    add_heatmap = TRUE,
    verbose = TRUE
)
Sys.time() - t0 #

# layout(cbind(1,2))
# plotAverage(E003_psa_no0) # cool
# plotAverage(E003_psa_w0)

myCols <- list(
    c("white", "red", "black"),
    c("white", "blue", "black"),
    c("white", "darkgreen", "black"),
    c("white", "red", "black"),
    c("white", "magenta", "black")
)
#
# plotHeatmap(E003_psa_no0, clstmethod = "none", sortrows = "increasing",
#             include = c(F, F, T, F, F), colvec = myCols)
#
# pdf(file = "../../../plots/E003_seqplots_heatmap.pdf", width = 11.69, height = 8.27)
# plotHeatmap(E003_psa_w0, clstmethod = "none", sortrows = "increasing",
#             include = c(F, F, T, F, F), colvec = myCols, raster = TRUE,
#             pointsize = 18, ln.v = FALSE)
# dev.off()
#
# pdf(file = "../../../plots/E003_seqplots_heatmap_no0.pdf", width = 11.69, height = 8.27)
# plotHeatmap(E003_psa_no0, clstmethod = "none", sortrows = "increasing",
#             include = c(F, F, T, F, F), colvec = myCols, raster = TRUE,
#             pointsize = 18, ln.v = FALSE)
# dev.off()

library(dplyr)
library(readr)

roadmap <- list(
    pc = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.pc"),
    nc = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.nc"),
    rb = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.rb")
)

lapply(roadmap, dim)
lapply(roadmap, head)
roadmap$all <- do.call(rbind, roadmap)

lapply(roadmap, function(x) anyDuplicated(x$gene_id))

E003_geneExp <- select(roadmap$all, gene_id, E003)

extractValuesFormPsa <- function(myPsa) {
    myPsaL <- myPsa$data[[1]] # one bed file only
    res <- lapply(myPsaL, function(x) x$heatmap)
    geneNames <- myPsa$annotations[[1]]$name # one bed file only
    geneNames <- sub("\\.[0-9]*", "", geneNames) # ENSG00000277232.2 to ENSG00000277232, should do nothing if already formated
    return(list(geneNames = geneNames, heatmaps = res))
}

matData <- extractValuesFormPsa(E003_psa_w0)

colnames(E003_geneExp) <- c("gene_id", "exp")

filterMatDataFromSortedExpData <- function(myMatData, myExpData) {
    colnames(myMatData$heatmaps[[1]]) <- paste0("FracMeth",   seq_len(ncol(myMatData$heatmaps[[1]])))
    colnames(myMatData$heatmaps[[2]]) <- paste0("CpG_sites",  seq_len(ncol(myMatData$heatmaps[[2]])))
    colnames(myMatData$heatmaps[[3]]) <- paste0("mCpG_dens",  seq_len(ncol(myMatData$heatmaps[[3]])))
    colnames(myMatData$heatmaps[[4]]) <- paste0("mCpG_ratio", seq_len(ncol(myMatData$heatmaps[[4]])))
    colnames(myMatData$heatmaps[[5]]) <- paste0("ReadCov",    seq_len(ncol(myMatData$heatmaps[[5]])))

    myMatDataT <- data.frame(
        gene_id = myMatData$geneNames,
        do.call(cbind, myMatData$heatmaps),
        stringsAsFactors = FALSE
    ) %>% tbl_df

    combTable <- inner_join(myExpData, myMatDataT, by = "gene_id")
    message(paste(
        "genes in expression table:    ", nrow(myExpData),
        "\ngenes in tracks data:         ", nrow(myMatDataT),
        "\ngenes common in both datasets:", nrow(combTable)
    ))
    combTable <- arrange(combTable, desc(exp))
    return(combTable)
}

combTable <- filterMatDataFromSortedExpData(matData, E003_geneExp)

refTable <- read_tsv(
    "../../../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.bed",
    col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)

addGeneTypeInfo <- function(myCombTable, myRefTable) {
    myCombTable <- inner_join(myRefTable, myCombTable, by = "gene_id")
    myCombTable <- arrange(myCombTable, desc(exp))
    myCombTable <- select(myCombTable, -chr, -start, -end, -score, -strand, -symbol)
    return(myCombTable)
}

combTable <- addGeneTypeInfo(combTable, refTable)

md <- combTable

library(plotrix)

plotDispersion <- function(myMat, bin = 5, colorPalette = rainbow, alphaForSe = 0.25) {
    # myMat <- combTable[, grep("CpG_sites", colnames(combTable))] %>% as.matrix
    myMats <- lapply(seq_len(bin), function(x) subset(myMat, ntile(seq_len(nrow(myMat)), bin) == x))
    myMeans <- lapply(myMats, function(x) colMeans(x, na.rm = TRUE))
    mySds <- lapply(myMats, function(x) apply(x, 2, function(y) sd(y, na.rm = TRUE)))
    mySes <- lapply(seq_along(myMeans), function(i) ifelse(myMeans[[i]] == 0, 0, mySds[[i]]/sqrt(length(myMeans[[i]]))))

    ymax <- max(
        sapply(
            seq_along(myMeans), function(i) {
                max(myMeans[[i]] + mySes[[i]])
            }
        )
    )
    xind <- seq_len(ncol(myMat))
    myColorPalette <- colorPalette(bin)

    plot(NA, xlim = range(xind), ylim = c(0, ymax), axes = FALSE, xlab = NA, ylab = NA) # \o/
    axis(1, at = c(xind[1], xind[(length(xind)+1)/2], xind[length(xind)]), labels = c("-2.5kb", "TSS", "+2.5kb"))
    axis(2)
    for(i in rev(seq_along(myMats))) {
        dispersion(xind, myMeans[[i]], mySes[[i]], type = "l", fill = adjustcolor(myColorPalette[i], alphaForSe))
        lines(xind, myMeans[[i]], type = "l", col = myColorPalette[i], lwd = 2)
    }
}


plotCombTable <- function(md,
                          colSuffix = c("FracMeth", "CpG_sites", "mCpG_dens", "mCpG_ratio", "ReadCov"),
                          colours = list(
                              c("white", "red", "black"),
                              c("white", "blue", "black"),
                              c("white", "green", "black"),
                              c("white", "red", "black"),
                              c("white", "magenta", "black")
                          ),
                          nbin = 5,
                          binPalette = colorRampPalette(c("magenta", "black", "green")),
                          raster = FALSE
                          ) {

    layout(matrix(c(1, 3, 4, 2, 3, 4, 5, 3, 3, 6:20), nrow = 3), widths = c(0.8, 0.15, 0.15, 1, 1, 1, 1, 1), heights = c(1, 0.12, 0.25))

    # expression plot
    par(mar = c(0, 2, 4, 0.5), mgp = c(1.4, 0.5, 0))
    plot(x = log10(md$exp +1), y = rev(seq_len(nrow(md))), type = "l", lwd = 3,
         axes = FALSE, xlab = NA, ylab = NA, ylim = c(1, nrow(md)),
         xaxs="i", yaxs="i")
    box()
    axis(3)
    axis(2, at = c(1, nrow(md)))
    mtext(side = 3, "log(FPKM+1)", line = 1.4, cex = 0.8)

    # geneType plot
    par(mar = c(0, 0.5, 4, 0))
    geneType <- rep(0, nrow(md))
    geneType[which(md$gene_type == "protein_coding")] <- 1
    geneType[which(md$gene_type == "processed_pseudogene")] <- 2
    image(t(rev(geneType)), col = c("blue", "green", "grey"), axes = FALSE, breaks = seq(-0.5, 2.5, by = 1))
    box()
    mtext(side = 3, "gene\ntype", line = 0.5,  cex = 0.8, las = 2)
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("top", legend = c("protein coding", "other", "pseudogene"), fill = c ("green", "blue", "grey"), bty = "n")

    # expression boxplot
    myExprs <- lapply(seq_len(nbin), function(x) log10(subset(md$exp, ntile(seq_len(nrow(md)), nbin) == x) + 1))
    par(mar = c(2.5, 2.5, 0, 1.5))
    boxplot(myExprs, col = binPalette(nbin), ylab = "log(FPKM+1)", xlab = "bin", pch = 19)

    # bin side legend
    par(mar = c(0, 0.5, 4, 0))
    image(t(rev(ntile(seq_len(nrow(md)), nbin))), col = binPalette(nbin), axes = FALSE)
    box()
    mtext(side = 3, "bin", line = 0.5,  cex = 0.8, las = 2)

    # heatmaps, col legend, mean plot
    for (i in seq_along(colSuffix)) {
        data1 <- md[, grep(colSuffix[i],colnames(md))] %>% as.matrix
        dataForAveragePlot <- data1

        # http://stackoverflow.com/questions/20977477/how-to-assign-a-specific-color-to-na-in-an-image-plot
        col <- colorRampPalette(colours[[i]])(255)
        zlim <- c(0, 1)
        if(i == 2) zlim <- c(0, 20)
        if(i == 3) zlim <- c(1, 4)
        if(i == 5) zlim <- c(0, 30)
        saved_zlim <- zlim
        newz.na <- zlim[2]+(zlim[2]-zlim[1])/length(col) # new z for NA
        newz.outside <- zlim[2]+2*(zlim[2]-zlim[1])/length(col) # new z for values outside zlim
        data1[which(is.na(data1>zlim[2]))] <- newz.na # we affect newz.outside
        data1[data1 > zlim[2]] <- newz.outside # same for newz.na
        zlim[2] <- zlim[2]+2*(zlim[2]-zlim[1])/length(col) # we finally extend the z limits to include the two new values
        col <- c(col, "grey", "black") # we construct the new color range by including: na.color and outside.color

        # heatmap
        par(mar = c(0, 1, 4, 1.5))
        image(
            t(data1[rev(seq_len(nrow(data1))), ]),
            col = col,
            zlim = zlim,
            useRaster = raster,
            axes = FALSE
        )
        axis(1, at = c(0, 0.5, 1), labels = c("-2.5kb", "TSS", "+2.5kb"))
        mtext(side = 3, colSuffix[i], line = 0.5,  cex = 0.8)

        # legend
        par(mar = c(2.2, 1, 2.2, 1.5))
        image(
            matrix(seq(saved_zlim[1], saved_zlim[2], length.out = 100)),
            col = colorRampPalette(colours[[i]])(255),
            useRaster = raster,
            axes = FALSE
        )
        box()
        axis(1, at = c(0, 1), labels = saved_zlim)

        # average plot
        par(mar = c(2, 1, 0, 1.5))
        plotDispersion(dataForAveragePlot, bin = nbin, colorPalette = binPalette, alphaForSe = 0.25)
    }

    layout(1)
}


matData <- extractValuesFormPsa(E003_psa_w0)
combTable <- filterMatDataFromSortedExpData(matData, E003_geneExp)
combTable <- addGeneTypeInfo(combTable, refTable)
pdf(file = "../../../plots/E003_my_heatmap.pdf", width = 11.69, height = 8.27)
plotCombTable(combTable, raster = TRUE)
dev.off()

pc <- filter(combTable, gene_type == "protein_coding")
pg <- filter(combTable, gene_type == "processed_pseudogene")
ot <- filter(combTable, !gene_type %in% c("protein_coding", "processed_pseudogene"))

pdf(file = "../../../plots/E003_my_heatmap_byGeneType.pdf", width = 11.69, height = 8.27)
plotCombTable(pc, raster = TRUE)
plotCombTable(ot, raster = TRUE)
plotCombTable(pg, raster = TRUE)
dev.off()

matData <- extractValuesFormPsa(E003_psa_no0)
combTable <- filterMatDataFromSortedExpData(matData, E003_geneExp)
combTable <- addGeneTypeInfo(combTable, refTable)
pdf(file = "../../../plots/E003_my_heatmap_no0.pdf", width = 11.69, height = 8.27)
plotCombTable(combTable, raster = TRUE)
dev.off()


