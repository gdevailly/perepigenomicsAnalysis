library(here)
library(tidyverse)
library(furrr); plan(multiprocess(workers = availableCores() - 2))
library(parallel)

# rename files
metadata <- read_tsv("perepigenomics/data/availableByFeature.tsv")
cell_dic <- read_tsv("metadata_33samples.txt")


# WGBS
wgbs_trs <- set_names(cell_dic$id, cell_dic$short)

fe <- new.env()
fe$rename_WGBS <- function(filename, trs) {
    load(here::here("perepigenomics/data/Rdata/", filename))
    byFeatureData <- map(byFeatureData, function(x) { # 30 s
        dplyr::mutate(x, cell_type = trs[cell_type])
    })
    save(byFeatureData, file = here::here("perepigenomics/data/Rdata/", filename))
    rm(byFeatureData)
    return(paste0(filename, " is done."))
}

t0 <- Sys.time()
future_map_chr(
    filter(metadata, assay == "WGBS")$file,
    ~fe$rename_WGBS(.x, trs = wgbs_trs)
)
Sys.time() - t0


# Dnase is already ok
filter(metadata, assay == "Dnase")

# HisMods
filter(metadata, !(assay %in% c("WGBS", "Dnase")))
fe <- new.env()
fe$rename_HisMod <- function(filename) {
    load(here("perepigenomics/data/Rdata/", filename))
    byFeatureData <- map(byFeatureData, function(x) { # 30 s
        mutate(x, cell_type = substr(cell_type, 1, 4))
    })
    save(byFeatureData, file = here("perepigenomics/data/Rdata/", filename))
    rm(byFeatureData)
    return(paste0(filename, " is done."))
}

t0 <- Sys.time()
res <- future_map_chr(
    filter(metadata, !(assay %in% c("WGBS", "Dnase")))$file,
    fe$rename_HisMod
)
Sys.time() - t0


# fe <- new.env()
# fe$rename <- function(filename) {
#     x <- load(here("perepigenomics/data/Rdata/", filename))
#     if (x == "tmp") {
#         byFeatureData <- get(x)
#         save(byFeatureData, file = here("perepigenomics/data/Rdata/", filename))
#     }
#     return(paste0(filename, " is done."))
# }
#
# t0 <- Sys.time()
# res <- future_map_chr(
#     metadata$file,
#     fe$rename
# )
# Sys.time() - t0


