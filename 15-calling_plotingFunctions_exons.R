setwd("~/mnt/genotoul_grp/guillaume/cascade")

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")

library(dplyr)
library(plotrix)
library(seqplots)
library(readr)
library(svglite)

metadata <- read_tsv("data/wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

epms <- read_rds("~/work/projects/cascade/data/Rdata/innerExonEpm.rds")

metadata$id[!metadata$id %in% colnames(epms)] # missing E008, E017, E021, E022
metadata <- filter(metadata, id %in% colnames(epms))

epms <- as.data.frame(epms)
epms$exon_location <- rownames(epms)

exonInfoPath <- "~/work/projects/cascade/data/inner_exon_epms.bed"
preffix <- "~/mnt/genotoul/work/projects/cascade/"

# exon TPM plot ------------------

# png version --------------------
t0 <- Sys.time() # 11 minutes
for(i in seq_len(nrow(metadata))){
    dataForPlot <- extractAndPrepareDataForExons(
        metadata$id[i],
        exonInfoPath,
        epms,
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
        file = paste0(preffix, "appPlots/exons/wgbs/byTpm/", metadata$id[i], ".png"),
        width = 8.3, height = 5.8, units = "in", res = 300, pointsize = 13
    )
    plotMetricAndProfile(
        dataForPlot,
        nbin = 5,
        zlims  = zlims,
        raster = TRUE,
        title = "",
        withGeneType = FALSE,
        expTransf = function(x) log10(x + 1),
        xlabels = c("-1kb", "Exon start", "+1kb"),
        expLim = c(0, 4),
        exp.text = "log10(TPM + 1)",
        reversedZOrder = TRUE,
        tints  = c("red"        , "blue"      , "purple"      , "sienna4"   )
    )
    dev.off()

    message(paste(metadata$id[i], "done!"))
}
Sys.time() - t0
