library(tidyverse)
library(furrr); plan(multiprocess(workers = availableCores() - 2))
library(parallel)

metadata <- read_tsv("perepigenomics/data/availableByFeature.tsv")

loadData <- function(md, what = "TSS", path = "perepigenomics/data/Rdata/") {
    files <- filter(md, feature == what)$file %>%
        set_names(filter(md, feature == what)$assay)
    ret <- furrr::future_map(seq_along(files), function(x) {
        load(here(path, files[x]))
        byFeatureData
    })
    set_names(ret, names(files))
}

t0 <- Sys.time() # 17s
TSS <- loadData(metadata, what = "TSS")
Sys.time() - t0

map(TSS, ~paste(colnames(.x[[1]]), collapse = " ")) %>% table(deparse.level = 0)
map_int(TSS, length) %>% table()

table_data <- mclapply(
    seq_along(TSS[[1]]),
    function(i) {
        gene_data <- map(TSS, i)
        gene_data[[1]] <- select(gene_data[[1]], cell_type, exp, mCpG_ratio)
        for (j in seq(2, length(gene_data))) {
            tmp <- select(gene_data[[j]], cell_type, HisMod)
            colnames(tmp)[2] <- names(gene_data)[j]
            gene_data[[j]] <- tmp
        }
        reduce(gene_data, full_join, by = "cell_type")
    },
    mc.cores = 12
)

names(table_data) <- names(TSS[[1]])

save(table_data, file = here("data/full_model_tss.RData"))

# New session --------------------
library(tidyverse)
library(furrr); plan(multiprocess(workers = availableCores() - 2))
library(parallel)

load("data/full_model_tss.RData")

layout(rbind(1, 2))
as.matrix(select(table_data[[1]], -cell_type)) %>%
    apply(1, function(x) sum(is.na(x))/length(x)) %>%
    sort() %>% plot(type = "l", ylim = c(0, 1), ylab = "Fraction of NA", xlab = "Cell types")
as.matrix(select(table_data[[1]], -cell_type)) %>%
    apply(2, function(x) sum(is.na(x))/length(x)) %>%
    sort() %>% plot(type = "l", ylim = c(0, 1), ylab = "Fraction of NA", xlab = "Epigenetic mark")


lm(exp ~ 1 + mCpG_ratio + Dnase, data = table_data[[1]])
lm(exp ~ 1 + mCpG_ratio + Dnase, data = table_data[[1]], na.action = na.exclude)

# keep only full column -------------
filtdata <- select(table_data[[1]], exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)
naniar::vis_miss(filtdata)
myformula <- paste("exp ~ 1 +", paste(colnames(filtdata)[-1], collapse = " + "))

fe <- new.env()
fe$get_stats <- function(x) {
    filtdata <- select(x, exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)
    filtdata$exp <- log10(x$exp + 1) # log10 +1 transformation of expression data
    default_output <- tibble(
        term = c("intercept", "mCpG_ratio", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"),
        statistic = NA_real_,
        p.value = NA_real_
    )
    if(any(is.nan(filtdata$mCpG_ratio))) { # NaN = region with no CpG
        return(default_output)
    }
    myformula <- paste("exp ~ 1 +", paste(colnames(filtdata)[-1], collapse = " + "))
    mlm <- lm(as.formula(myformula), data = filtdata, na.action = na.exclude) %>%
        broom::tidy() %>%
        mutate(term = sub("(Intercept)", "intercept", term, fixed = TRUE))
    if (nrow(mlm) != 7 || ncol(mlm) != 5) { # too many 0s, etc. on some genes, no coverage / missassembly
        return(default_output)
    }
    select(mlm, term, estimate, p.value)
}

stat_data <- future_map_dfr(table_data, fe$get_stats)
stat_data <- mutate(stat_data, gene = rep(names(table_data), each = 7)) %>% select(gene, everything())

p <- stat_data %>%
    filter(term != "intercept") %>%
    mutate(slope = case_when(
        p.value <= 0.01 & estimate > 0 ~ "Positive",
        p.value <= 0.01 & estimate < 0 ~ "Negative",
        p.value >  0.01                 ~ "N.S."      ,
        TRUE                            ~ "N.D."
    )) %>%
    filter(slope %in% c("Positive", "Negative")) %>%
    mutate(term = sub("mCpG_ratio", "DNAme", term, fixed = TRUE)) %>%
    mutate(term = factor(term, levels = c("DNAme", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))) %>%
    ggplot(aes(x = term, fill = slope)) +
    geom_bar(position = position_dodge()) +
    labs(x = "Mark", y = "Number of genes", fill = "Slope:", title = "Linear regression between gene expression\nand epigenetic mark at TSS") +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "cividis") +
    coord_cartesian(ylim = c(0, 2500)) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 20, hjust = 0.8))
ggsave(p, filename = "plots/global_lm_6marks_tss.png", width = 7, height = 5)

genemd <- read_tsv("perepigenomics/data/gene_list.tsv")
long <- filter(genemd, length_type == "long")$ensg

p <- stat_data %>%
    filter(gene %in% long) %>%
    filter(term != "intercept") %>% # TODO : ask Yann how to interpret intercept
    mutate(slope = case_when(
        p.value <= 0.01 & estimate > 0 ~ "Positive",
        p.value <= 0.01 & estimate < 0 ~ "Negative",
        p.value >  0.01                 ~ "N.S."      ,
        TRUE                            ~ "N.D."
    )) %>%
    filter(slope %in% c("Positive", "Negative")) %>%
    mutate(term = sub("mCpG_ratio", "DNAme", term, fixed = TRUE)) %>%
    mutate(term = factor(term, levels = c("DNAme", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))) %>%
    ggplot(aes(x = term, fill = slope)) +
    geom_bar(position = position_dodge()) +
    labs(x = "Mark", y = "Number of genes", fill = "Slope:", title = "Linear regression between long gene expression\nand epigenetic mark at TSS") +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "cividis") +
    coord_cartesian(ylim = c(0, 1000)) +
    scale_y_continuous(sec.axis = sec_axis(~ . / length(long), labels = scales::percent)) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 20, hjust = 0.8))
ggsave(p, filename = "plots/global_lm_6marks_tss_long_only.png", width = 7, height = 5)


# TTS -----------------
library(tidyverse)
library(furrr); plan(multiprocess(workers = availableCores() - 2))
library(parallel)

metadata <- read_tsv("perepigenomics/data/availableByFeature.tsv")

loadData <- function(md, what = "TSS", path = "perepigenomics/data/Rdata/") {
    files <- filter(md, feature == what)$file %>%
        set_names(filter(md, feature == what)$assay)
    ret <- furrr::future_map(seq_along(files), function(x) {

        load(paste0(path, files[x]))
        byFeatureData
    })
    set_names(ret, names(files))
}

t0 <- Sys.time() # 17s
TTS <- loadData(metadata, what = "TTS")
Sys.time() - t0

map(TTS, ~paste(colnames(.x[[1]]), collapse = " ")) %>% table(deparse.level = 0)
map_int(TTS, length) %>% table()

table_data <- mclapply(
    seq_along(TTS[[1]]),
    function(i) {
        gene_data <- map(TTS, i)
        gene_data[[1]] <- select(gene_data[[1]], cell_type, exp, mCpG_ratio)
        for (j in seq(2, length(gene_data))) {
            tmp <- select(gene_data[[j]], cell_type, HisMod)
            colnames(tmp)[2] <- names(gene_data)[j]
            gene_data[[j]] <- tmp
        }
        reduce(gene_data, full_join, by = "cell_type")
    },
    mc.cores = 12
)

names(table_data) <- names(TTS[[1]])

naniar::vis_miss(table_data[[1]])
filtdata <- select(table_data[[1]], exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)
naniar::vis_miss(filtdata)

fe <- new.env()
fe$get_stats <- function(x) {
    filtdata <- select(x, exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)
    filtdata$exp <- log10(x$exp + 1) # log10 +1 transformation of expression data
    default_output <- tibble(
        term = c("intercept", "mCpG_ratio", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"),
        statistic = NA_real_,
        p.value = NA_real_
    )
    if(any(is.nan(filtdata$mCpG_ratio))) { # NaN = region with no CpG
        return(default_output)
    }
    myformula <- paste("exp ~ 1 +", paste(colnames(filtdata)[-1], collapse = " + "))
    mlm <- lm(as.formula(myformula), data = filtdata, na.action = na.exclude) %>%
        broom::tidy() %>%
        mutate(term = sub("(Intercept)", "intercept", term, fixed = TRUE))
    if (nrow(mlm) != 7 || ncol(mlm) != 5) { # too many 0s, etc. on some genes, no coverage / missassembly
        return(default_output)
    }
    select(mlm, term, statistic, p.value)
}

stat_data <- future_map_dfr(table_data, fe$get_stats)
stat_data <- mutate(stat_data, gene = rep(names(table_data), each = 7)) %>% select(gene, everything())

p <- stat_data %>%
    filter(term != "intercept") %>%
    mutate(slope = case_when(
        p.value <= 0.01 & statistic > 0 ~ "Positive",
        p.value <= 0.01 & statistic < 0 ~ "Negative",
        p.value >  0.01                 ~ "N.S."      ,
        TRUE                            ~ "N.D."
    )) %>%
    filter(slope %in% c("Positive", "Negative")) %>%
    mutate(term = sub("mCpG_ratio", "DNAme", term, fixed = TRUE)) %>%
    mutate(term = factor(term, levels = c("DNAme", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))) %>%
    ggplot(aes(x = term, fill = slope)) +
    geom_bar(position = position_dodge()) +
    labs(x = "Mark", y = "Number of genes", fill = "Slope:", title = "Linear regression between gene expression\nand epigenetic mark at TTS") +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "cividis") +
    coord_cartesian(ylim = c(0, 2500)) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 20, hjust = 0.8))
ggsave(p, filename = "plots/global_lm_6marks_tts.png", width = 7, height = 5)

genemd <- read_tsv("perepigenomics/data/gene_list.tsv")
long <- filter(genemd, length_type == "long")$ensg

p <- stat_data %>%
    filter(gene %in% long) %>%
    filter(term != "intercept") %>%
    mutate(slope = case_when(
        p.value <= 0.01 & statistic > 0 ~ "Positive",
        p.value <= 0.01 & statistic < 0 ~ "Negative",
        p.value >  0.01                 ~ "N.S."      ,
        TRUE                            ~ "N.D."
    )) %>%
    filter(slope %in% c("Positive", "Negative")) %>%
    mutate(term = sub("mCpG_ratio", "DNAme", term, fixed = TRUE)) %>%
    mutate(term = factor(term, levels = c("DNAme", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))) %>%
    ggplot(aes(x = term, fill = slope)) +
    geom_bar(position = position_dodge()) +
    labs(x = "Mark", y = "Number of genes", fill = "Slope:", title = "Linear regression between long gene expression\nand epigenetic mark at TTS") +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "cividis") +
    coord_cartesian(ylim = c(0, 5000)) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 20, hjust = 0.8))
ggsave(p, filename = "plots/global_lm_6marks_tts_long_only.png", width = 7, height = 5)

# exon Psi ---------------------
library(tidyverse)
library(furrr); plan(multiprocess(workers = availableCores() - 2))
library(parallel)

metadata <- read_tsv("perepigenomics/data/availableByFeature.tsv")

loadData <- function(md, what = "TSS", path = "perepigenomics/data/Rdata/") {
    files <- dplyr::filter(md, feature == what)$file %>%
        set_names(dplyr::filter(md, feature == what)$assay)
    ret <- furrr::future_map(seq_along(files), function(x) {
        message(x)
        load(paste0(path, files[x]))
        byFeatureData
    })
    set_names(ret, names(files))
}

t0 <- Sys.time() # 17s
exon <- loadData(metadata, what = "exonPsi")
Sys.time() - t0

map(exon, ~paste(colnames(.x[[1]]), collapse = " ")) %>% table(deparse.level = 0)
map_int(exon, length) %>% table()

table_data <- mclapply(
    seq_along(exon[[1]]),
    function(i) {
        message(i)
        gene_data <- map(exon, i)
        if(names(gene_data)[1] != "WGBS") {
            gene_data <- c(gene_data[length(gene_data)], gene_data[seq(1, length(gene_data) - 1)])
            # gene_data[[1]]$cell_type <- c("E003", "E004", "E005", "E006", "E007", "E011", "E012", "E013", "E016", "E050", "E065", "E066", "E079", "E084", "E085", "E094", "E095", "E096", "E097", "E098", "E100", "E104", "E105", "E106", "E109", "E112", "E113")
        }
        gene_data[[1]] <- select(gene_data[[1]], cell_type, exp, mCpG_ratio)
        for (j in seq(2, length(gene_data))) {
            tmp <- select(gene_data[[j]], cell_type, HisMod)
            colnames(tmp)[2] <- names(gene_data)[j]
            gene_data[[j]] <- tmp
        }
        reduce(gene_data, full_join, by = "cell_type")
    },
    mc.cores = 12
)

names(table_data) <- names(exon[[1]])

naniar::vis_miss(table_data[[1]])
filtdata <- select(table_data[[1]], exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)
naniar::vis_miss(filtdata)

fe <- new.env()
fe$get_stats <- function(x) {
    filtdata <- select(x, exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)

    default_output <- tibble(
        term = c("intercept", "mCpG_ratio", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"),
        statistic = NA_real_,
        p.value = NA_real_
    )
    if(any(is.nan(filtdata$mCpG_ratio))) { # NaN = region with no CpG
        return(default_output)
    }
    myformula <- paste("exp ~ 1 +", paste(colnames(filtdata)[-1], collapse = " + "))
    mlm <- lm(as.formula(myformula), data = filtdata, na.action = na.exclude) %>%
        broom::tidy() %>%
        mutate(term = sub("(Intercept)", "intercept", term, fixed = TRUE))
    if (nrow(mlm) != 7 || ncol(mlm) != 5) { # too many 0s, etc. on some genes, no coverage / missassembly
        return(default_output)
    }
    select(mlm, term, statistic, p.value)
}

stat_data <- future_map_dfr(table_data, fe$get_stats)
stat_data <- mutate(stat_data, gene = rep(names(table_data), each = 7)) %>% select(gene, everything())

p <- stat_data %>%
    filter(term != "intercept") %>%
    mutate(slope = case_when(
        p.value <= 0.01 & statistic > 0 ~ "Positive",
        p.value <= 0.01 & statistic < 0 ~ "Negative",
        p.value >  0.01                 ~ "N.S."      ,
        TRUE                            ~ "N.D."
    )) %>%
    filter(slope %in% c("Positive", "Negative")) %>%
    mutate(term = sub("mCpG_ratio", "DNAme", term, fixed = TRUE)) %>%
    mutate(term = factor(term, levels = c("DNAme", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))) %>%
    ggplot(aes(x = term, fill = slope)) +
    geom_bar(position = position_dodge()) +
    labs(x = "Mark", y = "Number of exons", fill = "Slope:", title = "Linear regression between middle exon Psi\nand epigenetic mark at middle exon") +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "cividis") +
    coord_cartesian(ylim = c(0, 1000)) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 20, hjust = 0.8))
ggsave(p, filename = "plots/global_lm_6marks_exons_Psi.png", width = 7, height = 5)


# exon Tpm ---------------------
library(tidyverse)
library(furrr); plan(multiprocess(workers = availableCores() - 2))
library(parallel)

metadata <- read_tsv("perepigenomics/data/availableByFeature.tsv")

loadData <- function(md, what = "TSS", path = "perepigenomics/data/Rdata/") {
    files <- dplyr::filter(md, feature == what)$file %>%
        set_names(dplyr::filter(md, feature == what)$assay)
    ret <- furrr::future_map(seq_along(files), function(x) {
        message(x)
        load(paste0(path, files[x]))
        byFeatureData
    })
    set_names(ret, names(files))
}

t0 <- Sys.time() # 8s
exon <- loadData(metadata, what = "exonTpm")
Sys.time() - t0

map(exon, ~paste(colnames(.x[[1]]), collapse = " ")) %>% table(deparse.level = 0)
map_int(exon, length) %>% table()

table_data <- mclapply(
    seq_along(exon[[1]]),
    function(i) {
        gene_data <- map(exon, i)
        if(names(gene_data)[1] != "WGBS") {
            gene_data <- c(gene_data[length(gene_data)], gene_data[seq(1, length(gene_data) - 1)])
            # gene_data[[1]]$cell_type <- c("E003", "E004", "E005", "E006", "E007", "E011", "E012", "E013", "E016", "E050", "E065", "E066", "E079", "E084", "E085", "E094", "E095", "E096", "E097", "E098", "E100", "E104", "E105", "E106", "E109", "E112", "E113")
        }
        gene_data[[1]] <- select(gene_data[[1]], cell_type, exp, mCpG_ratio)
        for (j in seq(2, length(gene_data))) {
            tmp <- select(gene_data[[j]], cell_type, HisMod)
            colnames(tmp)[2] <- names(gene_data)[j]
            gene_data[[j]] <- tmp
        }
        reduce(gene_data, full_join, by = "cell_type")
    },
    mc.cores = 12
)

names(table_data) <- names(exon[[1]])

naniar::vis_miss(table_data[[1]])
filtdata <- select(table_data[[1]], exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)
naniar::vis_miss(filtdata)

fe <- new.env()
fe$get_stats <- function(x) {
    filtdata <- select(x, exp, mCpG_ratio, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)
    filtdata$exp <- log10(x$exp + 1) # log10 +1 transformation of expression data
    # ^ comment if Psi

    default_output <- tibble(
        term = c("intercept", "mCpG_ratio", "H3K27me3", "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3"),
        statistic = NA_real_,
        p.value = NA_real_
    )
    if(any(is.nan(filtdata$mCpG_ratio))) { # NaN = region with no CpG
        return(default_output)
    }
    myformula <- paste("exp ~ 1 +", paste(colnames(filtdata)[-1], collapse = " + "))
    mlm <- lm(as.formula(myformula), data = filtdata, na.action = na.exclude) %>%
        broom::tidy() %>%
        mutate(term = sub("(Intercept)", "intercept", term, fixed = TRUE))
    if (nrow(mlm) != 7 || ncol(mlm) != 5) { # too many 0s, etc. on some genes, no coverage / missassembly
        return(default_output)
    }
    select(mlm, term, statistic, p.value)
}

stat_data <- future_map_dfr(table_data, fe$get_stats)
stat_data <- mutate(stat_data, gene = rep(names(table_data), each = 7)) %>% select(gene, everything())

p <- stat_data %>%
    filter(term != "intercept") %>%
    mutate(slope = case_when(
        p.value <= 0.01 & statistic > 0 ~ "Positive",
        p.value <= 0.01 & statistic < 0 ~ "Negative",
        p.value >  0.01                 ~ "N.S."      ,
        TRUE                            ~ "N.D."
    )) %>%
    filter(slope %in% c("Positive", "Negative")) %>%
    mutate(term = sub("mCpG_ratio", "DNAme", term, fixed = TRUE)) %>%
    mutate(term = factor(term, levels = c("DNAme", "H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3"))) %>%
    ggplot(aes(x = term, fill = slope)) +
    geom_bar(position = position_dodge()) +
    labs(x = "Mark", y = "Number of exons", fill = "Slope:", title = "Linear regression between middle exon TPM\nand epigenetic mark at middle exon") +
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option = "cividis") +
    coord_cartesian(ylim = c(0, 1000)) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 20, hjust = 0.8))
ggsave(p, filename = "plots/global_lm_6marks_exons_Tpm.png", width = 7, height = 5)
