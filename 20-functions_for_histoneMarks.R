library(Repitools)
library(dplyr)
library(GenomicRanges)


# readAlign file are kind of bed files found on roadmap .ftp
importReadAlignAsGranges <- function(
    filePath,
    chromosomes_to_keep = paste0("chr", 1:22),
    col_types = "ciicic"
) {
    read_tsv(paste0("histone/", filePath), col_names = FALSE, col_types = col_types, progress = FALSE) %>% # import is long, file is huge
        dplyr::select(-X4, -X5) %>%
        dplyr::rename(chr = X1, start = X2, end = X3, strand = X6) %>%
        dplyr::filter(chr %in% chromosomes_to_keep) %>%
        with(., GRanges(chr, IRanges(start, end), strand = strand)) %>%
        return()
}

# take a GRanges of reads, an annotation table (see featureScores documentation) and return a normalised matrix (fpkm values)
extractFeatureScoreMatrix <- function(
    genoRange,
    annoTable,
    up = 2500,
    down = 2500,
    freq = 100,
    endFeature = FALSE,
    ...
) {
    fs <- featureScores(
        genoRange,
        annoTable,
        up = up,
        down = down,
        freq = freq,
        use.strand = FALSE, # reads must not be on same strand as feature
        ...
    ) %>%
        tables %>%
        .[[1]]
    fs <- fs * 1E6 * 1000/freq # transform score into fpkm values
    if(endFeature) {
        fs <- fs[, rev(seq_len(ncol(fs)))]
    }
    return(fs)
}


# prepare DNAse and control data for plot
# warnings are for absent mapability in Repitools::featureScores, can be ignored
prepareDNAseData <- function(myCellCode, tssAnno, metadataTable, expressionTable, endFeature = FALSE, ...) {

    signalMat <- importReadAlignAsGranges(
        dplyr::filter(metadataTable, cellCode == myCellCode)$DNAse_file
    ) %>% extractFeatureScoreMatrix(tssAnno, verbose = FALSE, endFeature = endFeature, ...)
    controlMat <- importReadAlignAsGranges(
        dplyr::filter(metadataTable, cellCode == myCellCode)$Control_file
    ) %>% extractFeatureScoreMatrix(tssAnno, verbose = FALSE, endFeature = endFeature, ...)

    signalDf <- data.frame(
        gene_id = rownames(signalMat),
        signalMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(signalDf) <- c("gene_id", paste0("DNAse", colnames(signalDf)[-1]))
    controlDf <- data.frame(
        gene_id = rownames(controlMat),
        controlMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(controlDf) <- c("gene_id", paste0("Control", colnames(controlDf)[-1]))

    return(
        inner_join(
            select_(expressionTable, "gene_id", myCellCode),
            signalDf,
            by = "gene_id"
        ) %>%
            inner_join(
                controlDf,
                by = "gene_id"
            ) %>%
            rename_(exp = myCellCode) %>%
            arrange(desc(exp))
    )
}

prepareDNAseDataExons <- function(myCellCode, annoTable, metadataTable, metricTable, endFeature = FALSE, ...) {

    signalMat <- importReadAlignAsGranges(
        dplyr::filter(metadataTable, cellCode == myCellCode)$DNAse_file
    ) %>% extractFeatureScoreMatrix(annoTable, verbose = FALSE, endFeature = endFeature, ...)
    controlMat <- importReadAlignAsGranges(
        dplyr::filter(metadataTable, cellCode == myCellCode)$Control_file
    ) %>% extractFeatureScoreMatrix(annoTable, verbose = FALSE, endFeature = endFeature, ...)

    signalDf <- data.frame(
        gene_id = rownames(signalMat),
        signalMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(signalDf) <- c("name", paste0("DNAse", colnames(signalDf)[-1]))
    controlDf <- data.frame(
        gene_id = rownames(controlMat),
        controlMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(controlDf) <- c("name", paste0("Control", colnames(controlDf)[-1]))

    return(
        inner_join(
            select_(metricTable, "name", myCellCode),
            signalDf,
            by = "name"
        ) %>%
            inner_join(
                controlDf,
                by = "name"
            ) %>%
            rename_(exp = myCellCode) %>%
            arrange(desc(exp))
    )
}

prepareHisModData <- function(i, tssAnno, metadataTable, expressionTable, endFeature = FALSE, ...) {

    signalMat <- importReadAlignAsGranges(
        metadataTable[i,]$file
    ) %>% extractFeatureScoreMatrix(tssAnno, verbose = FALSE, endFeature = endFeature, ...)
    controlMat <- importReadAlignAsGranges(
        metadataTable[i,]$input_file
    ) %>% extractFeatureScoreMatrix(tssAnno, verbose = FALSE, endFeature = endFeature, ...)

    signalDf <- data.frame(
        gene_id = rownames(signalMat),
        signalMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(signalDf) <- c("gene_id", paste0("HisMod", colnames(signalDf)[-1]))
    controlDf <- data.frame(
        gene_id = rownames(controlMat),
        controlMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(controlDf) <- c("gene_id", paste0("Input", colnames(controlDf)[-1]))

    return(
        inner_join(
            select_(expressionTable, "gene_id",  metadataTable[i,]$cellCode),
            signalDf,
            by = "gene_id"
        ) %>%
            inner_join(
                controlDf,
                by = "gene_id"
            ) %>%
            rename_(exp = metadataTable[i,]$cellCode) %>%
            arrange(desc(exp))
    )
}

prepareHisModDataExons <- function(i, annoTable, metadataTable, metricTable, endFeature = FALSE, ...) {

    signalMat <- importReadAlignAsGranges(
        metadataTable[i,]$file
    ) %>% extractFeatureScoreMatrix(annoTable, verbose = FALSE, endFeature = endFeature, ...)
    controlMat <- importReadAlignAsGranges(
        metadataTable[i,]$input_file
    ) %>% extractFeatureScoreMatrix(annoTable, verbose = FALSE, endFeature = endFeature, ...)

    signalDf <- data.frame(
        gene_id = rownames(signalMat),
        signalMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(signalDf) <- c("name", paste0("HisMod", colnames(signalDf)[-1]))
    controlDf <- data.frame(
        gene_id = rownames(controlMat),
        controlMat
    ) %>% as_data_frame %>%
        mutate(gene_id = sub("\\.[0-9]*", "", gene_id))
    colnames(controlDf) <- c("name", paste0("Input", colnames(controlDf)[-1]))

    return(
        inner_join(
            select_(metricTable, "name",  metadataTable[i,]$cellCode),
            signalDf,
            by = "name"
        ) %>%
            inner_join(
                controlDf,
                by = "name"
            ) %>%
            rename_(exp = metadataTable[i,]$cellCode) %>%
            arrange(desc(exp))
    )
}



