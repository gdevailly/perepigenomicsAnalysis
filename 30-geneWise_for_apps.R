setwd("/home/gdevailly/work/projects/cascade/perepigenomics/data/Rdata/")
library(tidyverse)

# rename files
files <- list.files() %>% grep(".RData$", ., value = TRUE)

# fast, but can easily be parallelized
t0 <- Sys.time() # 3 minutes
walk(files, function(file) {
    myVar <- load(file)
    byFeatureData <- get(myVar)
    rm(myVar)
    save(byFeatureData, file = file)
    message(paste(file, "has been renamed!"))
})
Sys.time() - t0


# metadata contruction

byFeatureMd <- tibble(
    file = list.files(),
    feature = case_when(
        grepl("exonFpkm" , file, fixed = TRUE) ~ "exonTpm",
        grepl("exonPsi" , file, fixed = TRUE) ~ "exonPsi",
        grepl("_tes_", file, fixed = TRUE) ~ "TTS",
        grepl("_tss_", file, fixed = TRUE) ~ "TSS",
        grepl("exonWiseData_cascade.RData", file, fixed = TRUE) ~ "exonTpm",
        TRUE ~ NA_character_
    ),
    assay = case_when(
        grepl("_cascade.RData", file, fixed = TRUE) ~ "WGBS",
        file == "geneWiseData_exonPsi_wgbs.RData" ~ "WGBS",
        TRUE ~ strsplit(file, "_", fixed = TRUE) %>%
            map_chr(last) %>%
            sub(".RData", "", .)
    )
)

table(byFeatureMd$feature) # ok, 24 of each
table(byFeatureMd$assay) # ok, 4 everywhere

write_tsv(byFeatureMd, "../availableByFeature.tsv")
