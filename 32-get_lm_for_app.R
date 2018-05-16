library(here)
library(tidyverse)

source(here("Rscripts", "11-geneWiseFunctions.R"))
source(here("Rscripts", "29-geneWiseFunctions_hisMods.R"))

byFeatureMd <- read_tsv(here("perepigenomics", "data", "availableByFeature.tsv"))

t0 <- Sys.time() # 2h
walk(
    seq_len(nrow(byFeatureMd)),
    function(i) {
        load(here("perepigenomics", "data", "Rdata", byFeatureMd$file[i]))

        if (byFeatureMd$assay[i] == "WGBS") {
            modelTable <- getLmAndSd(byFeatureData)
        } else {
            modelTable <- getLmAndSd_dnase(byFeatureData)
        }

        save(
            modelTable,
            file = here("perepigenomics", "data", "modelData", paste0("model_", byFeatureMd$file[i]))
        )

        message(paste("Done for", byFeatureMd$file[i], " !"))
    }
)
Sys.time() - t0
