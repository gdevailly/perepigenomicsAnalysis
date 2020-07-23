
extractGeneWiseDataForHistone <- function(
    geneName,
    dataList,
    windows = c("X.500", "X.400", "X.300", "X.200", "X.100", "X0", "X100", "X200", "X300", "X400", "X500") # middle 11 windows when windows go from 1 to 51
) {
    data_frame(
        cell_type = names(dataList),
        exp = vapply(
            dataList,
            function(x) filter(x, gene_id == geneName)$exp,
            0
        ),
        HisMod = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("HisMod", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        Input = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("Input", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        )
    )
}

extractGeneWiseDataForDnase <- function(
    geneName,
    dataList,
    windows = c("X.500", "X.400", "X.300", "X.200", "X.100", "X0", "X100", "X200", "X300", "X400", "X500") # middle 11 windows when windows go from 1 to 51
) {
    data_frame(
        cell_type = names(dataList),
        exp = vapply(
            dataList,
            function(x) filter(x, gene_id == geneName)$exp,
            0
        ),
        HisMod = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("DNAse", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        ),
        Input = vapply(
            dataList,
            function(x) {
                geneData <- filter(x, gene_id == geneName)
                geneData <- geneData[, colnames(geneData) %in% paste0("Control", windows)]
                rowMeans(geneData, na.rm = TRUE)
            },
            0
        )
    )
}


getLmAndSd_dnase <- function(listOfTbl, trans_func = function(x) log10(x + 1)) {
    data_frame(
        gene_id = names(listOfTbl),
        cv_exp = vapply(
            listOfTbl,
            function(x) sd(trans_func(x$exp), na.rm = TRUE)/mean(trans_func(x$exp), na.rm = TRUE),
            0
        ),
        cv_dnase = vapply(
            listOfTbl,
            function(x) sd(x$HisMod, na.rm = TRUE)/mean(x$HisMod, na.rm = TRUE),
            0
        ),
        sl_dnase = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$HisMod))) {
                    NA
                } else {
                    lm(trans_func(exp) ~ HisMod, x)$coefficients[2]
                }
            },
            0
        ),
        r2_dnase = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$HisMod))) {
                    NA
                } else {
                    summary(lm(trans_func(exp) ~ HisMod, x))$r.squared
                }
            },
            0
        ),
        cv_control = vapply(
            listOfTbl,
            function(x) sd(x$Input, na.rm = TRUE)/mean(x$Input, na.rm = TRUE),
            0
        ),
        sl_control = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$Input))) {
                    NA
                } else {
                    # lm(trans_func(x$exp) ~ Input, x)$coefficients[2]
                    lm(trans_func(exp) ~ Input, x)$coefficients[2]
                }
            },
            0
        ),
        r2_control = vapply(
            listOfTbl,
            function(x) {
                if (all(is.na(x$Input))) {
                    NA
                } else {
                    # summary(lm(trans_func(x$exp) ~ Input, x))$r.squared
                    summary(lm(trans_func(exp) ~ Input, x))$r.squared
                }
            },
            0
        )
    )
}
