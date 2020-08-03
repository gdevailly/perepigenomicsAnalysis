setwd("~/mnt/genotoul_grp/guillaume/cascade")

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/11-geneWiseFunctions.R")

library(dplyr)
library(plotrix)
library(seqplots)
library(readr)
library(svglite)
library(parallel)

metadata <- read_tsv("data/wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

psis <- readRDS("~/work/projects/cascade/data/Rdata/innerExonPsi.rds")
exon_location <- rownames(psis)
psis <- bind_cols(
    tibble(exon_location = exon_location),
    as_tibble(psis)
)

epms <- read_rds("~/work/projects/cascade/data/Rdata/innerExonEpm.rds")
epms <- as.data.frame(epms)
epms$exon_location <- rownames(epms)
epms$name <- epms$exon_location

metadata$id[!metadata$id %in% colnames(psis)] # missing "E008" "E017" "E021" "E022" "E024" "E053" "E054" "E058" "E070" "E071"
metadata <- filter(metadata, id %in% colnames(psis))

setwd("wgbs/roadmap")

exonInfoPath <- "~/work/projects/cascade/data/inner_exon_epms.bed"
exonTable <- read_tsv(exonInfoPath, col_names = FALSE, progress = FALSE)
colnames(exonTable) <- c("chr", "start", "end", "name", "score", "strand")

preffix <- "~/mnt/genotoul/work/projects/cascade/"


# psis ---------------
setwd("data/wgbs/roadmap")
exonInfoPath <- "~/work/projects/cascade/data/inner_exon_psis.bed"

t0 <- Sys.time()
dataForAllSamples <- mclapply(
    seq_len(nrow(metadata)),
    function(i) {
        message(i)
        dataForPlot <- extractAndPrepareDataForExons(
            metadata$id[i],
            exonInfoPath,
            psis,
            refgenome = "hg19",
            bin = 50L,
            rm0 = TRUE,
            xmin = 1000L, xmax = 1000L, type = "mf",
            add_heatmap = TRUE,
            verbose = FALSE
        )

        message(paste(metadata$id[i], "done!"))
        return(dataForPlot)
    },
    mc.cores = 6, mc.preschedule = FALSE
)
Sys.time() - t0 # 2 minutes
names(dataForAllSamples) <- metadata$id

t0 <- Sys.time()
exonWiseData <- mclapply(
    dataForAllSamples[[1]]$exon_location,
    function(x) extractExonWiseDataFor(x, dataForAllSamples, windows = 19:23), # 19:23 -> -100 +100 bp
    mc.cores = 14
)
Sys.time() - t0 # long, 25 minutes
names(exonWiseData) <- dataForAllSamples[[1]]$exon_location
byFeatureData <- exonWiseData
save(byFeatureData, file = "~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_exonPsi_wgbs.RData")


# EPMS------

t0 <- Sys.time()
dataForAllSamples <- mclapply(
    seq_len(nrow(metadata)),
    function(i) {
        dataForPlot <- extractAndPrepareDataForExons(
            metadata$id[i],
            exonInfoPath,
            epms,
            refgenome = "hg19",
            bin = 50L,
            rm0 = TRUE,
            xmin = 1000L, xmax = 1000L, type = "mf",
            add_heatmap = TRUE,
            verbose = FALSE
        )

        message(paste(metadata$id[i], "done!"))
        return(dataForPlot)
    },
    mc.cores = 6, mc.preschedule = FALSE
)
Sys.time() - t0 # 2 minutes
names(dataForAllSamples) <- metadata$id

t0 <- Sys.time()
exonWiseData <- mclapply(
    dataForAllSamples[[1]]$exon_location,
    function(x) extractExonWiseDataFor(x, dataForAllSamples, windows = 19:23), # 19:23 -> -100 +100 bp
    mc.cores = 14
)
Sys.time() - t0 # long, 25 minutes
names(exonWiseData) <- dataForAllSamples[[1]]$exon_location
byFeatureData <- exonWiseData
save(byFeatureData, file = "~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_exonTpm_wgbs.RData")


# exploratory plots ------------

load("../../Rdata/exonWiseData_cascade.RData")
library(cowplot)

n <- 1000
p <- plot_grid(
    exonWiseData[[n]] %>% ggplot(aes(x = CpG_sites   , y = log2(exp))) + geom_point(col = "darkred"  ) + geom_smooth(method = "lm"),
    exonWiseData[[n]] %>% ggplot(aes(x = mCpG_ratio  , y = log2(exp))) + geom_point(col = "darkblue" ) + geom_smooth(method = "lm"),
    exonWiseData[[n]] %>% ggplot(aes(x = mCpG_density, y = log2(exp))) + geom_point(col = "purple"   ) + geom_smooth(method = "lm"),
    exonWiseData[[n]] %>% ggplot(aes(x = coverage    , y = log2(exp))) + geom_point(col = "darkgreen") + geom_smooth(method = "lm"),
    ncol = 4
)
# save_plot("../../../plots/lm_exon1000.pdf", p, base_height = 4, base_width = 14)

# log2(exon_ratios) ---------------
t0 <- Sys.time()
modelTable <- getLmAndSd(exonWiseData, trans_func = log2)
Sys.time() - t0 # 9 minutes

save_plot("../../../plots/slopes_exon_inclusion_ratios.pdf", base_height = 4, base_width = 14,
          ggdraw() +
              draw_plot(myBoxplotFunc(modelTable, "CpG_sites"   , "red"      ), 0,    0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable, "mCpG_ratio"  , "blue"     ), 0.25, 0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable, "mCpG_density", "purple"   ), 0.5 , 0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable, "coverage"    , "darkgreen"), 0.75, 0, 0.25, 0.95) +
              draw_label("Middle exons, regression analysis", x = 0.5, y = 0.975)
)

# identity ---------------
t0 <- Sys.time()
modelTable_id <- getLmAndSd(exonWiseData, trans_func = function(x) return(x))
Sys.time() - t0 # 9 minutes

save_plot("../../../plots/slopes_exon_inclusion_ratios_identity.pdf", base_height = 4, base_width = 14,
          ggdraw() +
              draw_plot(myBoxplotFunc(modelTable_id, "CpG_sites"   , "red"      ), 0,    0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable_id, "mCpG_ratio"  , "blue"     ), 0.25, 0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable_id, "mCpG_density", "purple"   ), 0.5 , 0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable_id, "coverage"    , "darkgreen"), 0.75, 0, 0.25, 0.95) +
              draw_label("Middle exons, regression analysis", x = 0.5, y = 0.975)
)

pdf("../../../plots/lms_exon_inclusion_ratios_identifty.pdf", height = 4, width = 14)
myNs <- sample(seq_along(exonWiseData), size = 100)
for(n in myNs) {
    p <- plot_grid(
        exonWiseData[[n]] %>% ggplot(aes(x = CpG_sites   , y = exp)) + geom_point(col = "darkred"  ) + geom_smooth(method = "lm"),
        exonWiseData[[n]] %>% ggplot(aes(x = mCpG_ratio  , y = exp)) + geom_point(col = "darkblue" ) + geom_smooth(method = "lm"),
        exonWiseData[[n]] %>% ggplot(aes(x = mCpG_density, y = exp)) + geom_point(col = "purple"   ) + geom_smooth(method = "lm"),
        exonWiseData[[n]] %>% ggplot(aes(x = coverage    , y = exp)) + geom_point(col = "darkgreen") + geom_smooth(method = "lm"),
        ncol = 4
    )
    print(p)
}
dev.off()

# max = 1 ---------------
t0 <- Sys.time()
modelTable_max1 <- getLmAndSd(exonWiseData, trans_func = function(x) {
    x[which(x > 1)] <- 1
    return(x)
})
Sys.time() - t0 # 9 minutes

save_plot("../../../plots/slopes_exon_inclusion_ratios_max1.pdf", base_height = 4, base_width = 14,
          ggdraw() +
              draw_plot(myBoxplotFunc(modelTable_max1, "CpG_sites"   , "red"      ), 0,    0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable_max1, "mCpG_ratio"  , "blue"     ), 0.25, 0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable_max1, "mCpG_density", "purple"   ), 0.5 , 0, 0.25, 0.95) +
              draw_plot(myBoxplotFunc(modelTable_max1, "coverage"    , "darkgreen"), 0.75, 0, 0.25, 0.95) +
              draw_label("Middle exons, regression analysis", x = 0.5, y = 0.975)
)







