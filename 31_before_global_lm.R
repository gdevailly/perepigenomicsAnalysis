library(here)
library(tidyverse)
library(furrr); plan(multiprocess(workers = availableCores() - 2))
library(parallel)

# rename files
metadata <- read_tsv("perepigenomics/data/availableByFeature.tsv")
cell_dic <- read_tsv("metadata_33samples.txt")


# WGBS
wgbs_trs <- set_names(cell_dic$id, cell_dic$short)

cells27 <- c("E003", "E004", "E005", "E006", "E007", "E011", "E012", "E013",
            "E016", "E050", "E065", "E066", "E079", "E084", "E085", "E094",
            "E095", "E096", "E097", "E098", "E100", "E104", "E105", "E106",
            "E109", "E112", "E113")

# fe <- new.env()
# fe$rename_WGBS <- function(filename, trs) {
#     load(here::here("perepigenomics/data/Rdata/", filename))
#     byFeatureData <- map(byFeatureData, function(x) { # 30 s
#         dplyr::mutate(x, cell_type = trs[cell_type])
#     })
#     save(byFeatureData, file = here::here("perepigenomics/data/Rdata/", filename))
#     rm(byFeatureData)
#     return(paste0(filename, " is done."))
# }
#
# t0 <- Sys.time()
# future_map_chr(
#     filter(metadata, assay == "WGBS")$file,
#     ~fe$rename_WGBS(.x, trs = wgbs_trs)
# )
# Sys.time() - t0

load("perepigenomics/data/Rdata/geneWiseData_exonPsi_wgbs.RData", verbose = TRUE)
byFeatureData <- map(byFeatureData, function(x) {
             dplyr::mutate(x, cell_type = cells27)
        })
save(byFeatureData, file = "perepigenomics/data/Rdata/geneWiseData_exonPsi_wgbs.RData")
rm(byFeatureData)

load("perepigenomics/data/Rdata/geneWiseData_exonTpm_wgbs.RData", verbose = TRUE)
byFeatureData <- map(byFeatureData, function(x) {
    dplyr::mutate(x, cell_type = wgbs_trs[cell_type])
})
save(byFeatureData, file = "perepigenomics/data/Rdata/geneWiseData_exonTpm_wgbs.RData")
rm(byFeatureData)

load("perepigenomics/data/Rdata/geneWiseData_tes_cascade.RData", verbose = TRUE)
byFeatureData <- map(byFeatureData, function(x) {
    dplyr::mutate(x, cell_type = wgbs_trs[cell_type])
})
save(byFeatureData, file = "perepigenomics/data/Rdata/geneWiseData_tes_cascade.RData")
rm(byFeatureData)

load("perepigenomics/data/Rdata/geneWiseData_tss_cascade.RData", verbose = TRUE)
byFeatureData <- map(byFeatureData, function(x) {
    dplyr::mutate(x, cell_type = wgbs_trs[cell_type])
})
save(byFeatureData, file = "perepigenomics/data/Rdata/geneWiseData_tss_cascade.RData")
rm(byFeatureData)

# Dnase is already ok
filter(metadata, assay == "Dnase")
load("perepigenomics/data/Rdata/geneWiseData_tss_Dnase.RData", verbose = TRUE)

# HisMods
filter(metadata, !(assay %in% c("WGBS", "Dnase")))
fe <- new.env()
fe$rename_HisMod <- function(filename) {
    message(filename)
    load(here("perepigenomics/data/Rdata/", filename))
    byFeatureData <- byFeatureData[sapply(byFeatureData, function(x) !is.null(x))] # -> 17517 -> 16787 exons
    if(any(nchar(byFeatureData[[1]]$cell_type) > 4)) {
        byFeatureData <- map(byFeatureData, function(x) { # 30 s
            mutate(x, cell_type = substr(cell_type, 1, 4))
        })
        save(byFeatureData, file = here("perepigenomics/data/Rdata/", filename))
    }
    rm(byFeatureData)
    return(filename)
}

t0 <- Sys.time() # 8 mins
res <- future_map_chr(
    filter(metadata, !(assay %in% c("WGBS", "Dnase")))$file,
    fe$rename_HisMod
)
Sys.time() - t0

# geneWiseData_exonPsi_H2A.Z.RData
# geneWiseData_exonPsi_H2AK5ac.RData
# geneWiseData_exonPsi_H2BK120ac.RData
# Error in UseMethod("mutate_") :
#     pas de méthode pour 'mutate_' applicable pour un objet de classe "NULL"
# De plus : Warning message:
#     `mutate_()` is deprecated as of dplyr 0.7.0.
# Please use `mutate()` instead.
# See vignette('programming') for more help
# This warning is displayed once every 8 hours.
# Call `lifecycle::last_warnings()` to see where this warning was generated.
