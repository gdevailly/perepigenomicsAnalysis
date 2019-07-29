setwd("/groups2/joshi_grp/guillaume/cascade/")
setwd("~/mnt/joshi_grp2/cascade")

library(tidyverse)

availablePng <- list.files("perepigenomics/appPlots/", recursive = TRUE, full.names = FALSE)
availablePng <- strsplit(availablePng, "/", fixed = TRUE) %>%
    map_df(~tibble(
        feature = .x[1],
        assay = .x[2],
        gene_type = if_else(feature == "exons", NA_character_, .x[3]),
        metric = case_when(
            feature != "exons" ~ "gene FPKM",
            .x[3] == "byFpkm"  ~ "middle exons FPKM",
            .x[3] == "byPsi"   ~ "middle exons inclusion ratio",
            TRUE               ~ NA_character_
        ),
        cell = strsplit(.x[4], ".", fixed = TRUE)[[1]][1]
    ))

availablePng <- mutate(availablePng, feature_metric = case_when(
    feature == "tss" & metric == "gene FPKM" ~ "Transcription Start Sites (by gene FPKM)",
    feature == "tes" & metric == "gene FPKM" ~ "Transcription Termination Sites (by gene FPKM)",
    feature == "exons" & metric == "middle exons FPKM" ~ "Middle Exons (by exons FPKM)",
    feature == "exons" & metric == "middle exons inclusion ratio" ~ "Middle Exons (by inclusion ratio)",
    TRUE ~ NA_character_
))

map(availablePng, ~table(.x, useNA = "ifany"))
write_tsv(availablePng, "perepigenomics/data/availablePng.tsv")


myDictionary <- c(
    tss = "Transcription Start Sites",
    tes = "Transcription Termination Sites",
    dnase1 = "DNAse 1",
    wgbs = "WGBS",
    TSS = "Transcription Start Sites",
    TTS = "Transcription Termination Sites",
    exonFpkm = "Middle exons (by expression)",
    exonPsi = "Middle exons (by inclusion ratio)",
    short = "Short genes (<=1kb)",
    intermediate = "Genes between 1kb and 3kb long",
    long = "Long genes (> 3kb)",
    all = "All genes"
)
myDictionary <- tibble(inWord = names(myDictionary), outWord = myDictionary)
cells <- read_tsv("metadata_33samples.txt")
cell_dic <- tibble(
    inWord = cells$id,
    outWord = cells$name
)
myDictionary <- bind_rows(myDictionary, cell_dic)
write_tsv(myDictionary, "perepigenomics/data/dictionary.tsv")

