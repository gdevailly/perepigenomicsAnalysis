# setwd("~/work/projects/cascade/")
# setwd("~/work/project/perepigenomics_app")
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(purrr)

sym2genes <- readr::read_tsv("data/gene_list.tsv")
cell_code <- readr::read_tsv("data/dictionary.tsv")

tss <- list.files("data/Rdata/", pattern = "_tss_", full.names = TRUE)
names(tss) <- stringr::str_extract(tss, "[:alnum:]+\\.") %>%
    stringr::str_replace("\\.", "")

all_tables <- lapply(tss, function(x) {
    load(x, verbose = TRUE)
    byFeatureData
})
names(all_tables)[1] <- "WGBS"

tsg <- c("GCG", "INS", "MYOD1", "NPPB", "CLDN18")
tsgcode <- filter(sym2genes, symbol %in% tsg)$ensg

marks <- c("WGBS", "Dnase", "H3K4me3", "H3K27me3", "H3K27ac")

dfp <- map_dfr(marks, function(mark) {
    temp <- bind_rows(
        all_tables[[mark]][tsgcode],
        .id = "ensg"
    ) %>% mutate(mark = mark)
    if (mark == "WGBS") {
        select(temp, ensg, cell_type, mark, exp, mark_value = mCpG_ratio)
    } else {
        select(temp, ensg, cell_type, mark, exp, mark_value = HisMod)
    }
})

dfp <- left_join(dfp, sym2genes, by = "ensg") %>%
    left_join(cell_code, by = c("cell_type" = "inWord"))

dfp <- mutate(dfp,
    cell = if_else(
        outWord %in% c("Pancreas", "Lung", "Gastric", "Psoas_Muscle", "Right_Ventricle", "Left_Ventricle"),
        outWord,
        "Other"
    )
)

palette <- c("grey", rainbow(6, 0.5, 0.75))
names(palette) = unique(dfp$cell)

p <- ggplot(dfp, aes(x = mark_value, y = log10(exp + 1), colour = cell)) +
    geom_point() +
    facet_wrap(~ symbol + mark, scales = "free_x") +
    scale_color_manual(values = palette) + 
    labs(y = "log10(TPM + 1)", x = "Epigenetic mark", color = "Cell type:") +
    theme_bw()

ggsave("../perepigenomics_revision/tissue_specific.png", plot = p, width = 8, height = 8)
