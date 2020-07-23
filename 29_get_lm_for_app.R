library(here)
library(tidyverse)
library(furrr); plan(multiprocess(workers = 12))

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/11-geneWiseFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/29-geneWiseFunctions_hisMods.R")

byFeatureMd <- read_tsv(here("perepigenomics", "data", "availableByFeature.tsv"))

t0 <- Sys.time() # 53 minutes
future_map_chr(
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
        return(paste("Done for", byFeatureMd$file[i], " !"))
    }
)
Sys.time() - t0
