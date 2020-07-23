setwd("/media/gdevailly/SANS TITRE/inra/cascade")

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/11-geneWiseFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/20-functions_for_histoneMarks.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/29-geneWiseFunctions_hisMods.R")

library(readr)
library(parallel)
library(purrr)

# metadata building ----------------------
metadata <- read_tsv("wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

psis <- readRDS("~/work/projects/cascade/data/Rdata/innerExonPsi.rds")

exon_location <- rownames(psis)
psis <- bind_cols(
    data_frame(name = exon_location),
    as_data_frame(psis)
)

metadata$id[!metadata$id %in% colnames(psis)] # missing "E008" "E017" "E021" "E022" "E024" "E053" "E054" "E058" "E070" "E071"
metadata <- filter(metadata, id %in% colnames(psis))

exonInfoPath <- "~/work/projects/cascade/data/inner_exon_psis.bed"
exonTable <- read_tsv(exonInfoPath, col_names = FALSE, progress = FALSE)
colnames(exonTable) <- c("chr", "start", "end", "name", "score", "strand")

preffix <- "~/mnt/genotoul/work/projects/cascade/"


# histone ------------
hisFiles <- data_frame(
    file = list.files("histone/") %>%
        grep(".gz$", ., value = TRUE),
    name = gsub(".tagAlign.gz", "", file, fixed = TRUE)
)
hisFiles <- bind_cols(
    hisFiles,
    strsplit(hisFiles$name, "-") %>%
        map_df(~data_frame(cellCode = .x[1], ChIP = .x[2]))
)
his_md <- left_join(
    dplyr::filter(hisFiles, !ChIP %in% c("DNase", "DNaseControl", "Input")),
    dplyr::filter(hisFiles, ChIP == "Input") %>%
        dplyr::select(file, cellCode) %>%
        dplyr::rename(input_file = file),
    by = "cellCode"
)
his_md <- filter(his_md, cellCode %in% colnames(psis))

myHisMods <- names(table(his_md$ChIP)[table(his_md$ChIP) > 1])

# preparing data ------------------
# Histone ---------
preparDataFor <- function(thisHisMod) { # unpure

    mdForThisHisMod <- dplyr::filter(his_md, ChIP == thisHisMod)

    dataForThisHisMod <- lapply(
        seq_len(nrow(mdForThisHisMod)),
        function(i) {
            hisModData <- prepareHisModDataExons(
                i, annoTable = exonTable, metadataTable = mdForThisHisMod, metricTable = psis,
                up = 1000, down = 1000, freq = 50
            )
            message(paste(mdForThisHisMod$name[i], "done!"))
            return(hisModData)
        }
    )
    names(dataForThisHisMod) <- mdForThisHisMod$name
    dataForThisHisMod <- map(dataForThisHisMod, ~dplyr::rename(.x, gene_id = name))

    geneWiseData <- mclapply(
        dataForThisHisMod[[1]]$gene_id,
        function(x) extractGeneWiseDataForHistone(x, dataForThisHisMod, windows = c("X.100", "X.50", "X0", "X50", "X100")),
        mc.cores = 12
    )
    names(geneWiseData) <- dataForThisHisMod[[1]]$gene_id

    assign(thisHisMod, geneWiseData)
    save(list = thisHisMod, file = paste0(preffix, "Rdata/geneWiseData_exonPsi_", thisHisMod, ".RData"))

    return(NULL)
}

t0 <- Sys.time() # 11 h
walk(myHisMods, preparDataFor)
Sys.time() - t0


# DNAse1 ----------------------------
DNAse_md <- inner_join(
    dplyr::filter(hisFiles, ChIP  == "DNase") %>%
        dplyr::select(-name, -ChIP),
    dplyr::filter(hisFiles, ChIP  == "DNaseControl") %>%
        dplyr::select(-name, -ChIP),
    by = "cellCode"
) %>% dplyr::rename(DNAse_file = file.x, Control_file = file.y)

dataForDnase <- lapply(
    DNAse_md$cellCode,
    function(cellCode) {
        dnaseData <- prepareDNAseDataExons(cellCode, annoTable = exonTable, metadataTable = DNAse_md, metricTable = psis,
                                           up = 1000, down = 1000, freq = 50)
        message(paste(cellCode, "done!"))
        return(dnaseData)
    }
)
names(dataForDnase) <- DNAse_md$cellCode
dataForDnase <- map(dataForDnase, ~dplyr::rename(.x, gene_id = name))

t0 <- Sys.time()
Dnase <- mclapply(
    dataForDnase[[1]]$gene_id,
    function(x) extractGeneWiseDataForDnase(x, dataForDnase, windows = c("X.100", "X.50", "X0", "X50", "X100")),
    mc.cores = 12
)
Sys.time() - t0
names(Dnase) <- dataForDnase[[1]]$gene_id

save(Dnase, file = paste0(preffix, "Rdata/geneWiseData_exonPsi_Dnase.RData"))
