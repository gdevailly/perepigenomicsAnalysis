library(dplyr)
library(plotrix)
library(seqplots)
library(viridis)
library(raster)

# thanks https://stackoverflow.com/a/8189441
findMode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}


extractValuesFormPsa <- function(myPsa) {
    myPsaL <- myPsa$data[[1]] # one bed file only
    res <- lapply(myPsaL, function(x) x$heatmap)
    geneNames <- myPsa$annotations[[1]]$name # one bed file only
    geneNames <- sub("\\.[0-9]*", "", geneNames) # ENSG00000277232.2 to ENSG00000277232, should do nothing if already formated
    return(list(geneNames = geneNames, heatmaps = res))
}

# extract and format data using seqplots
# ... send to getPlotSetArray
# suggested:
# refgenome = "hg19",
# bin = 100L,
# rm0 = TRUE,
# xmin = 2500L, xmax = 2500L, type = "pf",
# add_heatmap = TRUE,
# verbose = TRUE
getPlotSetArrayFor <- function(rodmapCode, annotationPath, ...) {
    myBW <- paste0(rodmapCode, c(
        "_CpG_sites_hg19.bw",
        "_mCpG_ratio_hg19.bw",
        "_mCpG_density_hg19.bw",
        "_coverage_hg19.bw"
    ))
    myPsa <- getPlotSetArray(
        myBW,
        annotationPath,
        ...
    )
    myData <- extractValuesFormPsa(myPsa)
}

# extract signal from bigwig near features from path, and merge with gene expression fron table
# ... passed to seqplots::getPlotSetArray
extractAndPrepareDataFor <- function(roadmapCode, annotationPath, exprTable, verbose = TRUE, ...) {
    coverageData <- getPlotSetArrayFor(roadmapCode, annotationPath, verbose = verbose, ...)
    exprData <-  exprTable[, c("gene_id", roadmapCode)]
    colnames(exprData) <- c("gene_id", "exp")

    colnames(coverageData$heatmaps[[1]]) <- paste0("CpG_sites",     seq_len(ncol(coverageData$heatmaps[[1]])))
    colnames(coverageData$heatmaps[[2]]) <- paste0("mCpG_ratio",    seq_len(ncol(coverageData$heatmaps[[2]])))
    colnames(coverageData$heatmaps[[3]]) <- paste0("mCpG_density",  seq_len(ncol(coverageData$heatmaps[[3]])))
    colnames(coverageData$heatmaps[[4]]) <- paste0("coverage",      seq_len(ncol(coverageData$heatmaps[[4]])))

    myMatDataT <- data.frame(
        gene_id = coverageData$geneNames,
        do.call(cbind, coverageData$heatmaps),
        stringsAsFactors = FALSE
    ) %>% tbl_df

    combTable <- inner_join(exprData, myMatDataT, by = "gene_id")
    if (verbose) {
        message(paste(
            "genes in expression table:    ", nrow(exprData),
            "\ngenes in tracks data:         ", nrow(myMatDataT),
            "\ngenes common in both datasets:", nrow(combTable)
        ))
    }
    combTable <- arrange(combTable, desc(exp))
    return(combTable)
}

extractAndPrepareDataForWidth <- function(roadmapCode, annotationPath, ...) {
    coverageData <- getPlotSetArrayFor(roadmapCode, annotationPath, ...)

    colnames(coverageData$heatmaps[[1]]) <- paste0("CpG_sites",     seq_len(ncol(coverageData$heatmaps[[1]])))
    colnames(coverageData$heatmaps[[2]]) <- paste0("mCpG_ratio",    seq_len(ncol(coverageData$heatmaps[[2]])))
    colnames(coverageData$heatmaps[[3]]) <- paste0("mCpG_density",  seq_len(ncol(coverageData$heatmaps[[3]])))
    colnames(coverageData$heatmaps[[4]]) <- paste0("coverage",      seq_len(ncol(coverageData$heatmaps[[4]])))

    myMatDataT <- data.frame(
        gene_id = coverageData$geneNames,
        do.call(cbind, coverageData$heatmaps),
        stringsAsFactors = FALSE
    ) %>% tbl_df

    return(myMatDataT)
}

extractAndPrepareDataForExons <- function(roadmapCode, annotationPath, exonTable, verbose = TRUE, ...) {
    coverageData <- getPlotSetArrayFor(roadmapCode, annotationPath, verbose = verbose, ...)
    exonData <-  exonTable[, c("exon_location", roadmapCode)]
    colnames(exonData) <- c("exon_location", "exp")

    colnames(coverageData$heatmaps[[1]]) <- paste0("CpG_sites",     seq_len(ncol(coverageData$heatmaps[[1]])))
    colnames(coverageData$heatmaps[[2]]) <- paste0("mCpG_ratio",    seq_len(ncol(coverageData$heatmaps[[2]])))
    colnames(coverageData$heatmaps[[3]]) <- paste0("mCpG_density",  seq_len(ncol(coverageData$heatmaps[[3]])))
    colnames(coverageData$heatmaps[[4]]) <- paste0("coverage",      seq_len(ncol(coverageData$heatmaps[[4]])))

    myMatDataT <- data.frame(
        exon_location = coverageData$geneNames,
        do.call(cbind, coverageData$heatmaps),
        stringsAsFactors = FALSE
    ) %>% tbl_df

    combTable <- inner_join(exonData, myMatDataT, by = "exon_location")
    if (verbose) {
        message(paste(
            "genes in expression table:    ", nrow(exonData),
            "\ngenes in tracks data:         ", nrow(myMatDataT),
            "\ngenes common in both datasets:", nrow(combTable)
        ))
    }
    combTable <- arrange(combTable, desc(exp))
    return(combTable)
}

addGeneTypeInfo <- function(myCombTable, myRefTable) {
    myCombTable <- dplyr::inner_join(myRefTable, myCombTable, by = "gene_id")
    myCombTable <- dplyr::arrange(myCombTable, desc(exp))
    myCombTable <- dplyr::select(myCombTable, -chr, -start, -end, -score, -strand, -symbol)
    return(myCombTable)
}

# plot average profile plot, used by plotTraceAndKeyAndProfile ----------------------
plotDispersion <- function(
    myMat,
    bin = 5,
    bins = NULL, # bin is ignored when bins is pecified
    colorPalette = colorRampPalette(c("magenta", "black", "green")),
    alphaForSe = 0.25,
    xlabels = c("-2.5kb", "TSS", "+2.5kb"),
    reversedZOrder = FALSE
) {
    # myMat <- combTable[, grep("CpG_sites", colnames(combTable))] %>% as.matrix
    if(!is.null(bins)) {
        myMats <- lapply(unique(bins), function(x) subset(myMat, bins == x))
    } else {
        myMats <- lapply(seq_len(bin), function(x) subset(myMat, ntile(seq_len(nrow(myMat)), bin) == x))
    }
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
    if(!is.null(bins)) {
        myColorPalette <- colorPalette(length(unique(bins)))
    } else {
        myColorPalette <- colorPalette(bin)
    }

    plot(NA, xlim = range(xind), ylim = c(0, ymax), axes = FALSE, xlab = NA, ylab = NA) # \o/
    axis(1, at = c(xind[1], xind[(length(xind)+1)/2], xind[length(xind)]), labels = xlabels)
    axis(2)

    iter <- rev(seq_along(myMats))
    if (reversedZOrder) {
        iter <- rev(iter)
    }

    for(i in iter) {
        dispersion(xind, myMeans[[i]], mySes[[i]], type = "l", fill = adjustcolor(myColorPalette[i], alphaForSe))
        lines(xind, myMeans[[i]], type = "l", col = myColorPalette[i], lwd = 2)
    }
}

# make 3 plots: a trace heatmap plot, a colour key, and a profile plot
# layout to be call before the plot
# best stacked vertically
plotTraceAndKeyAndProfile <- function(
    myMat,
    tint = "red",
    naColour = "grey",
    blackEqualLow = FALSE,
    zlim = c(0, 1),
    raster = TRUE,
    xlabels = c("-2.5kb", "TSS", "+2.5kb"),
    title = "my data",
    profile.bin = 5,
    bins = NULL, # when specified, profile.bin is ignored
    profile.palette =  colorRampPalette(c("magenta", "black", "green")),
    profile.alphaForSe = 0.25,
    reversedZOrder = FALSE,
    npix_height = 250
) {
    if (tint == "viridis") {
        col <- rev(viridis(255))
        col <- c(col, naColour, col[255])
    } else {
        myColours <- c("white", tint, "black")
        if (blackEqualLow) myColours <- rev(myColours)
        col <- colorRampPalette(myColours)(255)
        col <- c(col, naColour, if(blackEqualLow) "white" else "black") # we construct the new color range by including: na.color and outside.color
    }

    # http://stackoverflow.com/questions/20977477/how-to-assign-a-specific-color-to-na-in-an-image-plot
    # https://stackoverflow.com/questions/11123152/function-for-resizing-matrices-in-r
    # matrix size reduction to avoid ploting artefact with image
    rr <- raster(myMat)
    extent(rr) <- extent(c(-180, 180, -90, 90))
    ss <- raster(nrow = npix_height, ncol = ncol(myMat))
    ss <- resample(rr, ss)

    saved_zlim <- zlim
    data1 <- as.matrix(ss)
    newz.na <- zlim[2]+(zlim[2]-zlim[1])/length(col) # new z for NA
    newz.outside <- zlim[2]+2*(zlim[2]-zlim[1])/length(col) # new z for values outside zlim
    data1[data1 > zlim[2]] <- newz.outside
    data1[which(is.na(data1))] <- newz.na
    zlim[2] <- zlim[2]+2*(zlim[2]-zlim[1])/length(col) # we finally extend the z limits to include the two new values

    oldMar <- par()$mar
    # heatmap
    par(mar = c(0, 1, 4, 1.5))

    image(
        t(data1[rev(seq_len(nrow(data1))), ]),
        col = col,
        zlim = zlim,
        useRaster = raster,
        axes = FALSE
    )
    axis(1, at = c(0, 0.5, 1), labels = xlabels)
    mtext(side = 3, title, line = 0.5,  cex = 0.8)

    # legend
    par(mar = c(2.2, 1, 2.2, 1.5))
    image(
        matrix(seq(saved_zlim[1], saved_zlim[2], length.out = 100)),
        col = if(tint == "viridis") rev(viridis(255)) else colorRampPalette(myColours)(255),
        useRaster = raster,
        axes = FALSE
    )
    box()
    axis(1, at = c(0, 1), labels = saved_zlim)

    # average plot
    par(mar = c(2, 1, 0, 1.5))
    if(!is.null(bins)) {
        plotDispersion(myMat, bins = bins, colorPalette = profile.palette, alphaForSe = profile.alphaForSe, xlabels = xlabels, reversedZOrder = reversedZOrder)
    } else {
        plotDispersion(myMat, bin = profile.bin, colorPalette = profile.palette, alphaForSe = profile.alphaForSe, xlabels = xlabels, reversedZOrder = reversedZOrder)
    }
    par(mar = oldMar)
}

# draw two plots, call layout first -----------------------
plotExpression <- function(
    myDF,
    profile.bin = 5,
    profile.palette =  colorRampPalette(c("magenta", "black", "green")),
    expLimLog10 = c(0, 4)
) {
    oldMar <- par()$mar

    # expression line plot
    par(mar = c(0, 2, 4, 0))
    plot(x = log10(myDF$exp +1), y = rev(seq_len(nrow(myDF))), type = "l", lwd = 3,
         axes = FALSE, xlab = NA, ylab = NA, ylim = c(1, nrow(myDF)), xlim = expLimLog10,
         xaxs="i", yaxs="i")
    box()
    axis(3)
    axis(2, at = c(1, nrow(myDF)))
    mtext(side = 3, "log10(FPKM+1)", line = 1.4, cex = 0.8)

    # expression boxplot
    myExprs <- lapply(seq_len(profile.bin), function(x) log10(subset(myDF$exp, ntile(seq_len(nrow(myDF)), profile.bin) == x) + 1))
    par(mar = c(2.5, 2.5, 0, 0))
    boxplot(myExprs, col = profile.palette(profile.bin), ylab = "log10(FPKM+1)", xlab = "bin", pch = 19, ylim = expLimLog10)
    par(mar = oldMar)
}

plotExpression2 <- function(
    myDF,
    profile.bin = 5,
    bins = NULL, # profile.bin is ignored when bins is specified
    profile.palette =  colorRampPalette(c("magenta", "black", "green")),
    expTransf = function(x) log10(x+1),
    expLim = c(0, 4),
    axis.text = "log10(FPKM+1)",
    abline.v = NULL
) {
    oldMar <- par()$mar

    # expression line plot
    par(mar = c(0, 2, 4, 0))
    plot(x = expTransf(myDF$exp), y = rev(seq_len(nrow(myDF))), type = "l", lwd = 3,
         axes = FALSE, xlab = NA, ylab = NA, ylim = c(1, nrow(myDF)), xlim = expLim,
         xaxs="i", yaxs="i")
    abline(v = abline.v, lty = 2, col = "grey")
    box()
    axis(3)
    axis(2, at = c(1, nrow(myDF)))
    mtext(side = 3, axis.text, line = 1.4, cex = 0.8)

    # expression boxplot
    par(mar = c(2.5, 2.5, 0, 0))
    if(!is.null(bins)) {
        myExprs <- lapply(unique(bins), function(x) expTransf(subset(myDF$exp, bins == x)))
        boxplot(myExprs, col = profile.palette(length(unique(bins))), ylab = axis.text, xlab = "bin", pch = 19, ylim = expLim)
    } else {
        myExprs <- lapply(seq_len(profile.bin), function(x) expTransf(subset(myDF$exp, ntile(seq_len(nrow(myDF)), profile.bin) == x)))
        boxplot(myExprs, col = profile.palette(profile.bin), ylab = axis.text, xlab = "bin", pch = 19, ylim = expLim)
    }
    par(mar = oldMar)

}

# draw two plots, call layout first -------------------
plotWidth <- function(
    myDF,
    profile.bin = 5,
    profile.palette =  colorRampPalette(c("magenta", "black", "green")),
    expLimLog10 = c(0, 4)
) {
    oldMar <- par()$mar

    # expression line plot
    par(mar = c(0, 2, 4, 0))
    plot(x = log10(myDF$end - myDF$start + 1), y = rev(seq_len(nrow(myDF))), type = "l", lwd = 3,
         axes = FALSE, xlab = NA, ylab = NA, ylim = c(1, nrow(myDF)), xlim = expLimLog10,
         xaxs="i", yaxs="i")
    box()
    axis(3)
    axis(2, at = c(1, nrow(myDF)))
    mtext(side = 3, "log10(width+1)", line = 1.4, cex = 0.8)

    # expression boxplot
    myWidth <- lapply(seq_len(profile.bin), function(x) log10(subset((myDF$end - myDF$start + 1), ntile(seq_len(nrow(myDF)), profile.bin) == x)))
    par(mar = c(2.5, 2.5, 0, 0))
    boxplot(myWidth, col = profile.palette(profile.bin), ylab = "log10(width+1)", xlab = "bin", pch = 19, ylim = expLimLog10)
    par(mar = oldMar)
}

# ---------------------
plotExpressionAndProfile <- function(
    myDF,
    nbin = 5,
    prefix = c("CpG_sites"  , "mCpG_ratio", "mCpG_density", "coverage"),
    labels = c("CpG density", "mCpG ratio", "mCpG density", "coverage"),
    tints  = c("red"        , "blue"      , "purple"      , "green"   ),
    zlims  = list(c(0, 20)  , c(0, 1)     , c(0, 4)       , c(0, 30)  ),
    naColour = "grey",
    geneTypePalette = rainbow,
    title = "",
    raster = TRUE,
    blackEqualLow = FALSE,
    expLimLog10 = c(0, 4),
    xlabels = c("-2.5kb", "TSS", "+2.5kb"),
    expressionLayoutWidth = c(0.8, 0.15, 0.15)
) {
    layout(
        matrix(c(1, 4, 2, 3, 4, 2, 5, 4, 6, 7:(6 + 3*length(prefix))), nrow = 3),
        widths = c(expressionLayoutWidth, rep(1, length(prefix))),
        heights = c(0.6, 0.15, 0.35)
    )
    oldMar <- par()$mar
    oldMgp <- par()$mgp

    par(mgp = c(1.4, 0.5, 0))
    # expression
    plotExpression(myDF, profile.bin = nbin, expLimLog10 = expLimLog10)

    # gene type
    par(mar = c(0, 0.8, 4, 0))
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
    image(t(rev(ntile(seq_len(nrow(myDF)), nbin))), col = colorRampPalette(c("magenta", "black", "green"))(nbin), axes = FALSE)
    axis(4, at = seq(0, 1, length.out = nbin + 1), labels = NA)
    box()
    mtext(side = 3, "bin", line = 0.5,  cex = 0.8, las = 2)
    plot.new()

    for (i in seq_along(prefix)) {
        myData <- myDF[, grep(prefix[i], colnames(myDF))] %>% as.matrix
        plotTraceAndKeyAndProfile(
            myData,
            tint = tints[i],
            blackEqualLow = blackEqualLow,
            zlim = zlims[[i]],
            raster = raster,
            title = labels[i],
            profile.bin = nbin,
            naColour = naColour,
            xlabels = xlabels
        )
    }

    mtext(title, line = -1, outer = TRUE)

    # reseting graphical parameters
    par(mar = oldMar, mgp = oldMgp)
    layout(1)
}
# ------------
plotWidthAndProfile <- function(
    myDF,
    nbin = 5,
    prefix = c("CpG_sites"  , "mCpG_ratio", "mCpG_density", "coverage"),
    labels = c("CpG density", "mCpG ratio", "mCpG density", "coverage"),
    tints  = c("red"        , "blue"      , "purple"      , "green"   ),
    zlims  = list(c(0, 20)  , c(0, 1)     , c(0, 4)       , c(0, 30)  ),
    naColour = "grey",
    geneTypePalette = rainbow,
    title = "",
    raster = TRUE,
    blackEqualLow = FALSE,
    expLimLog10 = c(0, 4),
    xlabels = c("-1kb", "exon center", "+1kb")
) {
    layout(
        matrix(c(1, 4, 2, 3, 4, 2, 5, 4, 6, 7:(6 + 3*length(prefix))), nrow = 3),
        widths = c(0.8, 0.15, 0.15, rep(1, length(prefix))),
        heights = c(0.6, 0.15, 0.35)
    )
    oldMar <- par()$mar
    oldMgp <- par()$mgp

    par(mgp = c(1.4, 0.5, 0))
    # expression
    plotWidth(myDF, profile.bin = nbin, expLimLog10 = expLimLog10)

    # gene type
    par(mar = c(0, 0.8, 4, 0))
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
    image(t(rev(ntile(seq_len(nrow(myDF)), nbin))), col = colorRampPalette(c("magenta", "black", "green"))(nbin), axes = FALSE)
    axis(4, at = seq(0, 1, length.out = nbin + 1), labels = NA)
    box()
    mtext(side = 3, "bin", line = 0.5,  cex = 0.8, las = 2)
    plot.new()

    for (i in seq_along(prefix)) {
        myData <- myDF[, grep(prefix[i], colnames(myDF))] %>% as.matrix
        plotTraceAndKeyAndProfile(
            myData,
            tint = tints[i],
            blackEqualLow = blackEqualLow,
            zlim = zlims[[i]],
            raster = raster,
            title = labels[i],
            profile.bin = nbin,
            naColour = naColour,
            xlabels = xlabels
        )
    }

    mtext(title, line = -2, outer = TRUE)

    # reseting graphical parameters
    par(mar = oldMar, mgp = oldMgp)
    layout(1)
}

# ---------------------------
plotMetricAndProfile <- function(
    myDF,
    nbin = 5,
    bins = NULL, # nbin is ignored if bin is specified. Bins should be a numeric vector of length nrow(myDF)
    prefix = c("CpG_sites"  , "mCpG_ratio", "mCpG_density", "coverage"),
    labels = c("CpG density", "mCpG ratio", "mCpG density", "coverage"),
    tints  = c("red"        , "blue"      , "purple"      , "green"   ),
    zlims  = list(c(0, 20)  , c(0, 1)     , c(0, 4)       , c(0, 30)  ),
    naColour = "grey",
    withGeneType = TRUE,
    geneTypePalette = rainbow,
    title = "",
    raster = TRUE,
    blackEqualLow = FALSE,
    expTransf = function(x) log10(x+1),
    expLim = c(0, 4),
    exp.text = "log10(FPKM+1)",
    xlabels = c("-2.5kb", "TSS", "+2.5kb"),
    reversedZOrder = FALSE,
    abline.v = NULL,
    expressionLayoutWidth = c(0.8, 0.15, 0.15),
    npix_height = 500
) {
    layout(
        matrix(c(1, 4, 2, 3, 4, 2, 5, 4, 6, 7:(6 + 3*length(prefix))), nrow = 3),
        widths = c(expressionLayoutWidth, rep(1, length(prefix))),
        heights = c(0.6, 0.15, 0.35)
    )
    oldMar <- par()$mar
    oldMgp <- par()$mgp

    par(mgp = c(1.4, 0.5, 0))
    # expression
    if(!is.null(bins)) {
        plotExpression2(myDF, bins = bins, expTransf = expTransf, expLim = expLim, axis.text = exp.text, abline.v = abline.v)
    } else {
        plotExpression2(myDF, profile.bin = nbin, expTransf = expTransf, expLim = expLim, axis.text = exp.text, abline.v = abline.v)
    }

    # gene type
    par(mar = c(0, 0.8, 4, 0))
    if(withGeneType) {
        geneType <- as.numeric(factor(myDF$gene_type))
        
        # compaction to avoid ploting artefacts
        geneTypeShort <- data_frame(gene_type = geneType, bin = ntile(seq_len(length(geneType)), npix_height)) %>%
            group_by(bin) %>%
            summarise(mode = findMode(gene_type)) %>%
            dplyr::select(mode) %>%
            unlist
        
        image(t(rev(geneTypeShort)), col = geneTypePalette(length(unique(myDF$gene_type))), axes = FALSE,
              breaks = seq(0.5, length(unique(myDF$gene_type)) + 0.5, by = 1))
        box()
        mtext(side = 3, "gene\ntype", line = 0.5,  cex = 0.8, las = 2)
        par(mar = c(0, 0, 0, 0))
        plot.new()
        ml <- legend("topright", legend = rep("", length(unique(myDF$gene_type))), fill = geneTypePalette(length(unique(myDF$gene_type))), bty = "n",
                     title = "gene type:", title.adj = -5)
        text(ml$text$x - 0.08, ml$text$y - 0.005, levels(factor(myDF$gene_type)), pos = 2)
    } else {
        plot.new()
        plot.new()
    }

    # bining
    par(mar = c(0, 0.8, 4, 0))
    if(!is.null(bins)) {
        
        # compaction to avoid ploting artefacts
        binShort <- data_frame(x = bins, bin = ntile(seq_len(length(bins)), npix_height)) %>%
            group_by(bin) %>%
            summarise(mode = findMode(x)) %>%
            dplyr::select(mode) %>%
            unlist
        
        image(t(rev(as.matrix(binShort))), col = colorRampPalette(c("magenta", "black", "green"))(length(unique(bins))), axes = FALSE)
        at <- c(1 - sapply(
            unique(bins),
            function(x) {
                which(bins == x)[1]
            }
        )/length(bins), 0)
        axis(4, at = at, labels = NA)
    } else {
        bins <- ntile(seq_len(nrow(myDF)), nbin)
        binShort <- data_frame(x = bins, bin = ntile(seq_len(length(bins)), npix_height)) %>%
            group_by(bin) %>%
            summarise(mode = findMode(x)) %>%
            dplyr::select(mode) %>%
            unlist
        
        image(t(rev(as.matrix(binShort))), col = colorRampPalette(c("magenta", "black", "green"))(nbin), axes = FALSE)
        axis(4, at = seq(0, 1, length.out = nbin + 1), labels = NA)
    }
    box()
    mtext(side = 3, "bin", line = 0.5,  cex = 0.8, las = 2)
    plot.new()

    for (i in seq_along(prefix)) {
        myData <- myDF[, grep(prefix[i], colnames(myDF))] %>% as.matrix
        if(!is.null(bins)) {
            plotTraceAndKeyAndProfile(
                myData,
                tint = tints[i],
                blackEqualLow = blackEqualLow,
                zlim = zlims[[i]],
                raster = raster,
                title = labels[i],
                bins = bins,
                naColour = naColour,
                xlabels = xlabels,
                reversedZOrder = reversedZOrder,
                npix_height = npix_height
            )
        } else {
            plotTraceAndKeyAndProfile(
                myData,
                tint = tints[i],
                blackEqualLow = blackEqualLow,
                zlim = zlims[[i]],
                raster = raster,
                title = labels[i],
                profile.bin = nbin,
                naColour = naColour,
                xlabels = xlabels,
                reversedZOrder = reversedZOrder,
                npix_height = npix_height
            )
        }
    }

    mtext(title, line = -1.3, outer = TRUE)

    # reseting graphical parameters
    par(mar = oldMar, mgp = oldMgp)
    layout(1)
}

# plot profiles only
plotProfileOnly <- function(
    myDF,
    nbin = 5,
    bins = NULL, # nbin is ignored if bin is specified. Bins should be a numeric vector of length nrow(myDF)
    prefix = c("CpG_sites"  , "mCpG_ratio", "mCpG_density", "coverage"),
    labels = c("CpG density", "mCpG ratio", "mCpG density", "coverage"),
    tints  = c("red"        , "blue"      , "purple"      , "green"   ),
    zlims  = list(c(0, 20)  , c(0, 1)     , c(0, 4)       , c(0, 30)  ),
    naColour = "grey",
    withGeneType = TRUE,
    geneTypePalette = rainbow,
    title = "",
    raster = TRUE,
    blackEqualLow = FALSE,
    xlabels = c("-2.5kb", "TSS", "+2.5kb"),
    reversedZOrder = FALSE,
    abline.v = NULL,
    npix_height = 500
) {
    layout(
        matrix(c(rep(1, 3), seq_len(3*length(prefix))+1), nrow = 3),
        widths = c(0.1, rep(1, length(prefix))),
        heights = c(0.65, 0.15, 0.30)
    )
    oldMar <- par()$mar
    oldMgp <- par()$mgp
    par(mar = c(0, 0, 0, 0))
    plot.new()
    for (i in seq_along(prefix)) {
        myData <- myDF[, grep(prefix[i], colnames(myDF))] %>% as.matrix
        if(!is.null(bins)) {
            plotTraceAndKeyAndProfile(
                myData,
                tint = tints[i],
                blackEqualLow = blackEqualLow,
                zlim = zlims[[i]],
                raster = raster,
                title = labels[i],
                bins = bins,
                naColour = naColour,
                xlabels = xlabels,
                reversedZOrder = reversedZOrder,
                npix_height = npix_height
            )
        } else {
            plotTraceAndKeyAndProfile(
                myData,
                tint = tints[i],
                blackEqualLow = blackEqualLow,
                zlim = zlims[[i]],
                raster = raster,
                title = labels[i],
                profile.bin = nbin,
                naColour = naColour,
                xlabels = xlabels,
                reversedZOrder = reversedZOrder,
                npix_height = npix_height
            )
        }
    }

    mtext(title, line = -1.3, outer = TRUE)

    # reseting graphical parameters
    par(mar = oldMar, mgp = oldMgp)
    layout(1)
}

