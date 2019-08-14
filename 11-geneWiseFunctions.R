library(cowplot)
library(broom)

extractGeneWiseDataFor <- function(
    geneName,
    dataList,
    windows = 21:31 # middle 11 windows when windows go from 1 to 51
) {
    data_frame(
        cell_type = names(dataList),
        exp = vapply(
            dataList,
            function(x) filter(x, gene_id == geneName)$exp,
            0
        ),
        CpG_sites = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("CpG_sites", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        mCpG_ratio = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("mCpG_ratio", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        mCpG_density = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("mCpG_density", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        coverage = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("coverage", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        )
    )
}

extractExonWiseDataFor <- function(
    geneName,
    dataList,
    windows = 17:25 # middle 11 windows when windows go from 1 to 41
) {
    data_frame(
        cell_type = names(dataList),
        exp = vapply(
            dataList,
            function(x) filter(x, exon_location == geneName)$exp,
            0
        ),
        CpG_sites = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, exon_location == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("CpG_sites", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        mCpG_ratio = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, exon_location == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("mCpG_ratio", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        mCpG_density = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, exon_location == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("mCpG_density", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        coverage = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, exon_location == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("coverage", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        )
    )
}

getLmAndSd <- function(listOfTbl, trans_func = function(x) log10(x + 1)) {
    data_frame(
        gene_id = names(listOfTbl),
        cv_exp = vapply(
            listOfTbl,
            function(x) sd(trans_func(x$exp), na.rm = TRUE)/mean(trans_func(x$exp), na.rm = TRUE),
            0
        ),
        cv_CpG_sites = vapply(
            listOfTbl,
            function(x) sd(x$CpG_sites, na.rm = TRUE)/mean(x$CpG_sites, na.rm = TRUE),
            0
        ),
        sl_CpG_sites = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$CpG_sites))) {
                    NA
                } else {
                    lm(trans_func(exp) ~ CpG_sites, x)$coefficients[2]
                }
            },
            0
        ),
        r2_CpG_sites = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$CpG_sites))) {
                    NA
                } else {
                    summary(lm(trans_func(exp) ~ CpG_sites, x))$r.squared
                }
            },
            0
        ),
        cv_mCpG_ratio = vapply(
            listOfTbl,
            function(x) sd(x$mCpG_ratio, na.rm = TRUE)/mean(x$mCpG_ratio, na.rm = TRUE),
            0
        ),
        sl_mCpG_ratio = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$mCpG_ratio))) {
                    NA
                } else {
                    lm(trans_func(exp) ~ mCpG_ratio, x)$coefficients[2]
                }
            },
            0
        ),
        r2_mCpG_ratio = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$mCpG_ratio))) {
                    NA
                } else {
                    summary(lm(trans_func(exp) ~ mCpG_ratio, x))$r.squared
                }
            },
            0
        ),
        cv_mCpG_density = vapply(
            listOfTbl,
            function(x) sd(x$mCpG_density, na.rm = TRUE)/mean(x$mCpG_density, na.rm = TRUE),
            0
        ),
        sl_mCpG_density = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$mCpG_density))) {
                    NA
                } else {
                    lm(trans_func(exp) ~ mCpG_density, x)$coefficients[2]
                }
            },
            0
        ),
        r2_mCpG_density = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$mCpG_density))) {
                    NA
                } else {
                    summary(lm(trans_func(exp) ~ mCpG_density, x))$r.squared
                }
            },
            0
        ),
        cv_coverage = vapply(
            listOfTbl,
            function(x) sd(x$coverage, na.rm = TRUE)/mean(x$coverage, na.rm = TRUE),
            0
        ),
        sl_coverage = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$coverage))) {
                    NA
                } else {
                    lm(trans_func(exp) ~ coverage, x)$coefficients[2]
                }
            },
            0
        ),
        r2_coverage = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$coverage))) {
                    NA
                } else {
                    summary(lm(trans_func(exp) ~ coverage, x))$r.squared
                }
            },
            0
        )
    )
}

cvFilter <- function(x, thresholds = c(0.25, 0.5, 0.75), values = c(0.125, 0.375, 0.625, 0.875)) {
    y <- rep(values[1], length(x))
    for(i in seq_along(thresholds)) {
        y[which(x > thresholds[i])] <- values[i + 1]
    }
    return(y)
}

myBoxplotFunc <- function(myTbl, myVarName, colour = "white") {
    newTbl <- mutate(myTbl, bin = cvFilter(myTbl[, paste0("r2_", myVarName)][[1]]))
    ylims <- sapply(unique(newTbl$bin), function(x) {
        newTbl2 <- filter(newTbl, bin == x)
        boxplot.stats(newTbl2[, paste0("sl_", myVarName)][[1]])$stats[c(1, 5)] * 1.05
    })
    ylims <- c(min(ylims[1,]), max(ylims[2,]))
    eff <- table(newTbl$bin)
    perc <- round(eff *100 / sum(eff), digits = 1)
    ann <- paste0(eff, "\n(", perc,"%)")
    return(ggplot(data = newTbl, aes_string("bin", paste0("sl_", myVarName), group = "bin")) + geom_boxplot(colour = colour, outlier.shape = NA) +
               coord_cartesian(ylim = ylims) +
               geom_hline(yintercept = 0) +
               labs(title = myVarName, x = bquote(R^2), y = "Slope") +
               annotate("text", x = unique(newTbl$bin), y = ylims[1], label = ann, vjust = 0))
}
