setwd("/media/gdevailly/DDCRCL/inra/cascade")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/20-functions_for_histoneMarks.R")

library(readr)
library(dplyr)
library(purrr)
library(tidyr)
library(parallel)

metadata <- read_tsv("~/mnt/genotoul_grp/guillaume/cascade/data/wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

psis <- readRDS("~/work/projects/cascade/data/Rdata/innerExonPsi.rds")

metadata$id[!metadata$id %in% colnames(psis)] # missing "E008" "E017" "E021" "E022" "E024" "E053" "E054" "E058" "E070" "E071"
metadata <- filter(metadata, id %in% colnames(psis))

psis <- as.data.frame(psis)
psis$exon_location <- rownames(psis)

exonInfoPath <- "~/work/projects/cascade/data/inner_exon_psis.bed"
preffix <- "~/mnt/genotoul/work/projects/cascade/"

# build histone marks table --------------
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

hisFiles <- filter(hisFiles, cellCode %in% colnames(psis))

his_md <- left_join(
    dplyr::filter(hisFiles, !ChIP %in% c("DNase", "DNaseControl", "Input")),
    dplyr::filter(hisFiles, ChIP == "Input") %>%
        dplyr::select(file, cellCode) %>%
        dplyr::rename(input_file = file),
    by = "cellCode"
)


# DNAseI -----------------------------
DNAse_md <- inner_join(
    dplyr::filter(hisFiles, ChIP  == "DNase") %>%
        dplyr::select(-name, -ChIP),
    dplyr::filter(hisFiles, ChIP  == "DNaseControl") %>%
        dplyr::select(-name, -ChIP),
    by = "cellCode"
) %>% dplyr::rename(DNAse_file = file.x, Control_file = file.y)


# exon Psi plot ------------------
exonTable <- read_tsv(exonInfoPath, col_names = FALSE, progress = FALSE)
colnames(exonTable) <- c("chr", "start", "end", "name", "score", "strand")
psis <- dplyr::rename(psis, name = exon_location)

# png versions -------------------
plotHisModDataForLine <- function(i) {
    hisModData <- prepareHisModDataExons(
        i, annoTable = exonTable, metadataTable = his_md, metricTable = psis,
        up = 1000, down = 1000, freq = 50
    )

    hisModName <- his_md[i, ]$ChIP
    cellID <-  his_md[i, ]$cellCode

    plotMetricAndProfile(
        hisModData,
        prefix = c("HisMod"  , "Input"),
        labels = c(hisModName, "Input"),
        tints  = c("sienna1" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        withGeneType = FALSE,
        expLim = c(0, 1.2),
        exp.text = "Exon PSI",
        expTransf = function(x) x,
        reversedZOrder = TRUE,
        abline.v = 1,
        xlabels = c("-1kb", "Exon start", "+1kb"),
        npix_height = 600
    )
    message(paste(his_md[i, ]$name, "is done!"))
}

t0 <- Sys.time() # 5.7h
lapply(
    137:242,
    #seq_len(nrow(his_md)),
    function(i) {
        hisMark <- his_md[i, ]$ChIP
        cellID  <- his_md[i, ]$cellCode
        if(!file.exists(paste0(preffix, "appPlots/exons/", hisMark, "/byPsi"))){
            dir.create(paste0(preffix, "appPlots/exons/", hisMark, "/byPsi"))
        }
        png(
            file = paste0(
                preffix, "appPlots/exons/", hisMark, "/byPsi/", cellID, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotHisModDataForLine(i)
        dev.off()
    } # big memory footprint
) %>% invisible()
Sys.time() - t0


# DNAseI -----------------------
plotDNAseDataForLine <- function(i) {
    dnaseData <- prepareDNAseDataExons(
        DNAse_md$cellCode[i], annoTable = exonTable, metadataTable = DNAse_md, metricTable = psis,
        up = 1000, down = 1000, freq = 50
    )

    cellID <- DNAse_md$cellCode[i]

    plotMetricAndProfile(
        dnaseData,
        prefix = c("DNAse"  , "Control"),
        labels = c("DNAse"  , "Control"),
        tints  = c("dodgerblue" , "grey"),
        zlims  = list(c(0, 2) , c(0, 2)),
        title = "",
        withGeneType = FALSE,
        expLim = c(0, 1.2),
        exp.text = "Exon PSI",
        expTransf = function(x) x,
        reversedZOrder = TRUE,
        abline.v = 1,
        xlabels = c("-1kb", "Exon start", "+1kb"),
        npix_height = 600
    )
    message(paste(cellID, "is done!"))
}

t0 <- Sys.time() # 31 minutes
lapply(
    seq_len(nrow(DNAse_md)),
    function(i) {
        cellID <- DNAse_md$cellCode[i]
        if(!file.exists(paste0(preffix, "appPlots/exons/dnase1/byPsi"))){
            dir.create(paste0(preffix, "appPlots/exons/dnase1"))
            dir.create(paste0(preffix, "appPlots/exons/dnase1/byPsi"))
        }
        png(
            file = paste0(
                preffix, "appPlots/exons/dnase1/byPsi/", cellID, ".png"
            ),
            width = 5.3, height = 5.8, units = "in", res = 300, pointsize = 13
        )
        plotDNAseDataForLine(i)
        dev.off()
    } # big memory footprint
)
Sys.time() - t0
