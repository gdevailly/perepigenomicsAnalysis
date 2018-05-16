setwd("/home/gdevailly/work/projects/cascade/perepigenomics/data/Rdata/")
library(tidyverse)

# rename files
files <- list.files() %>% grep(".RData$", ., value = TRUE)

# fast, but can easily be parallelized
walk(files, function(file) {
    myVar <- load(file)
    byFeatureData <- get(myVar)
    rm(myVar)
    save(byFeatureData, file = file)
    message(paste(file, "has been renamed!"))
})


# metadata contruction

byFeatureMd <- tibble(
    file = list.files(),
    feature = case_when(
        grepl("exon" , file, fixed = TRUE) ~ "exon",
        grepl("_tes_", file, fixed = TRUE) ~ "TES",
        grepl("_tss_", file, fixed = TRUE) ~ "TSS",
        TRUE ~ NA_character_
    ),
    assay = case_when(
        grepl("_cascade.RData", file, fixed = TRUE) ~ "WGBS",
        TRUE ~ strsplit(file, "_", fixed = TRUE) %>%
            map_chr(last) %>%
            sub(".RData", "", .)
    )
)

write_tsv(byFeatureMd, "../availableByFeature.tsv")
