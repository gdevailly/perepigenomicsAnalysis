setwd("/groups2/joshi_grp/guillaume/cascade/")
setwd("~/mnt/joshi_grp2/cascade")

library(tidyverse)

availablePng <- list.files("appPlots", recursive = TRUE, full.names = TRUE)
availablePng <- strsplit(availablePng, "/", fixed = TRUE) %>%
    map_df(~tibble(
        feature = .x[2],
        assay = .x[3],
        gene_type = if_else(feature == "exons", NA_character_, .x[4]),
        metric = case_when(
            feature != "exons" ~ "gene FPKM",
            .x[4] == "byFpkm"  ~ "exon FPKM",
            .x[4] == "byPsi"   ~ "exon inclusion ratio",
            TRUE               ~ NA_character_
        ),
        cell = strsplit(.x[5], ".", fixed = TRUE)[[1]][1]
    ))

availablePng <- mutate(availablePng, feature_metric = case_when(
    feature == "tss" & metric == "gene FPKM" ~ "TSS (gene FPKM)",
    feature == "tes" & metric == "gene FPKM" ~ "TES (gene FPKM)",
    feature == "exons" & metric == "exon FPKM" ~ "exons (exon FPKM)",
    feature == "exons" & metric == "exon inclusion ratio" ~ "exons (inclusion ratio)",
    TRUE ~ NA_character_
))

map(availablePng, ~table(.x, useNA = "ifany"))
write_tsv(availablePng, "data/availablePng.tsv")


myDictionary <- c(
    tss = "TSS",
    tes = "TES",
    dnase1 = "DNAse 1",
    wgbs = "WGBS"
)
myDictionary <- tibble(inWord = names(myDictionary), outWord = myDictionary)
cells <- read_tsv("metadata_33samples.txt")
cell_dic <- tibble(
    inWord = cells$id,
    outWord = cells$name
)
myDictionary <- bind_rows(myDictionary, cell_dic)
write_tsv(myDictionary, "perepigenomics/data/dictionary.tsv")

