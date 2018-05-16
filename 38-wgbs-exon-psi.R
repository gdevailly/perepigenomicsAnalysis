setwd("/groups2/joshi_grp/guillaume/cascade/data/")
source("../Rscripts/6-plotingFunctions.R")

library(dplyr)
library(plotrix)
library(seqplots)
library(readr)
library(svglite)

metadata <- read_tsv("wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

psis <- readRDS("Rdata/innerExonPsi.rds")
exon_location <- rownames(psis)
psis <- bind_cols(
    data_frame(exon_location = exon_location),
    as_data_frame(psis)
)

metadata$id[!metadata$id %in% colnames(psis)] # missing "E008" "E017" "E021" "E022" "E024" "E053" "E054" "E058" "E070" "E071"
metadata <- filter(metadata, id %in% colnames(psis))

setwd("wgbs/roadmap")
exonInfoPath <- "../../inner_exon_psis.bed"

# png plots ----------
t0 <- Sys.time() # 10 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataForExons(
        metadata$id[i],
        exonInfoPath,
        psis,
        refgenome = "hg19",
        bin = 50L,
        rm0 = TRUE,
        xmin = 1000L, xmax = 1000L, type = "pf",
        add_heatmap = TRUE,
        verbose = TRUE
    )

    zlim1 <- c(0, 20)
    zlim4 <- c(10, 50)
    if(metadata$id[i] %in% c(
        "E005", "E007", "E011", "E012", "E013", "E016", "E065", "E066", "E070", "E071", "E095", "E106"
    )) zlim4 <- c(30, 100)
    if(metadata$id[i] %in% c("E084", "E085")) {
        zlim1 <- c(0, 10)
        zlim4 <- c(0, 20)
    }
    zlims <- list(zlim1, c(0,1), c(1, 4), zlim4)

    png(
        file = paste0("../../../appPlots/exons/wgbs/byPsi/", metadata$id[i], ".png"),
        width = 8.3, height = 5.8, units = "in", res = 300, pointsize = 13
    )
    plotMetricAndProfile(
        dataForPlot,
        nbin = 5,
        zlims  = zlims,
        raster = TRUE,
        title = "",
        withGeneType = FALSE,
        expTransf = function(x) x,
        xlabels = c("-1kb", "Exon start", "+1kb"),
        expLim = c(0, 1.2),
        exp.text = "Exon PSI",
        reversedZOrder = TRUE,
        abline.v = 1,
        tints  = c("red"        , "blue"      , "purple"      , "sienna4"   )
    )
    dev.off()

    message(paste(metadata$id[i], "done!"))
}
Sys.time() - t0
