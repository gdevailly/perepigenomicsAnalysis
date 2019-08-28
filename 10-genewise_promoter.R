setwd("~/mnt/genotoul_grp/guillaume/cascade")

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/11-geneWiseFunctions.R")

library(readr)
library(parallel)

# metadata building
metadata <- read_tsv("data/wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

myProms <- "annotation/gencode.v29.annotation.hg19.middleTSStranscript.light.autosomes.bed"
salmonExp <- read_tsv("data/rnaseq/salmon_exp_genes_expressed.tsv")

refTable <- read_tsv(
  "annotation/gencode.v29.annotation.hg19.middleTSStranscript.bed",
  col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)
refTable <- filter(refTable, chr %in% paste0("chr", 1:22))
refTable <- filter(refTable, gene_id %in% salmonExp$gene_id)
metadata$id[!metadata$id %in% colnames(salmonExp)] # missing E008, E017, E021, E022, "E024" "E053" "E054" "E058" "E070" "E071"
metadata <- filter(metadata, id %in% colnames(salmonExp))

t0 <- Sys.time()
dataForAllSamples <- mclapply(
    seq_len(nrow(metadata)),
    function(i) {
        dataForPlot <- extractAndPrepareDataFor(
            metadata$id[i],
            myProms,
            salmonExp,
            refgenome = "hg19",
            bin = 100L,
            rm0 = TRUE,
            xmin = 2500L, xmax = 2500L, type = "pf",
            add_heatmap = TRUE,
            verbose = FALSE
        )
        dataForPlot <- addGeneTypeInfo(dataForPlot, refTable)
        message(paste(metadata$id[i], "done!"))
        return(dataForPlot)
    },
    mc.cores = 6, mc.preschedule = FALSE
)
Sys.time() - t0 # 2 minutes
names(dataForAllSamples) <- metadata$short

# gene extraction

t0 <- Sys.time()
geneWiseData <- mclapply(
    dataForAllSamples[[1]]$gene_id,
    function(x) extractGeneWiseDataFor(x, dataForAllSamples),
    mc.cores = 12
)
Sys.time() - t0 # long, ~30 minutes
names(geneWiseData) <- dataForAllSamples[[1]]$gene_id
save(geneWiseData, file = "~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_tss_cascade.RData")
load("~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_tss_cascade.RData")

n <- 33
p <- plot_grid(
    geneWiseData[[n]] %>% ggplot(aes(x = CpG_sites   , y = log10(exp + 1))) + geom_point(col = "darkred"  ) + geom_smooth(method = "lm"),
    geneWiseData[[n]] %>% ggplot(aes(x = mCpG_ratio  , y = log10(exp + 1))) + geom_point(col = "darkblue" ) + geom_smooth(method = "lm"),
    geneWiseData[[n]] %>% ggplot(aes(x = mCpG_density, y = log10(exp + 1))) + geom_point(col = "purple"   ) + geom_smooth(method = "lm"),
    geneWiseData[[n]] %>% ggplot(aes(x = coverage    , y = log10(exp + 1))) + geom_point(col = "darkgreen") + geom_smooth(method = "lm"),
    ncol = 4
)
save_plot("../../../plots/lm_ENSG00000152661.pdf", p, base_height = 4, base_width = 14)
getLmAndSd(geneWiseData[33])

t0 <- Sys.time()
modelTable <- getLmAndSd(geneWiseData) # default is log10(exp + 1)
Sys.time() - t0 # 12 minutes

cvFilter <- function(x, thresholds = c(0.25, 0.5, 0.75), values = c(0.125, 0.375, 0.625, 0.875)) {
    y <- rep(values[1], length(x))
    for(i in seq_along(thresholds)) {
        y[which(x > thresholds[i])] <- values[i + 1]
    }
    return(y)
}

myBoxplotFunc <- function(myTbl, myVarName, colour = "white") {
    newTbl <- mutate(myTbl, bin = cvFilter(myTbl[, paste0("r2_", myVarName)][[1]]))
    ylims <- sapply(unique(newTbl$bin), function(x) {
        newTbl2 <- filter(newTbl, bin == x)
        boxplot.stats(newTbl2[, paste0("sl_", myVarName)][[1]])$stats[c(1, 5)] * 1.05
    })
    ylims <- c(min(ylims[1,]), max(ylims[2,]))
    eff <- table(newTbl$bin)
    perc <- round(eff *100 / sum(eff), digits = 1)
    ann <- paste0(eff, "\n(", perc,"%)")
    return(ggplot(data = newTbl, aes_string("bin", paste0("sl_", myVarName), group = "bin")) + geom_boxplot(colour = colour, outlier.shape = NA) +
               coord_cartesian(ylim = ylims) +
               geom_hline(yintercept = 0) +
               labs(title = myVarName, x = bquote(R^2), y = "Slope") +
               annotate("text", x = unique(newTbl$bin), y = ylims[1], label = ann, vjust = 0))
}

myBoxplotFunc(modelTable, "mCpG_ratio", "blue")

save_plot("../../../plots/slopes_tss_allGenes.pdf", base_height = 4, base_width = 14,
    ggdraw() +
        draw_plot(myBoxplotFunc(modelTable, "CpG_sites"   , "red"      ), 0,    0, 0.25, 0.95) +
        draw_plot(myBoxplotFunc(modelTable, "mCpG_ratio"  , "blue"     ), 0.25, 0, 0.25, 0.95) +
        draw_plot(myBoxplotFunc(modelTable, "mCpG_density", "purple"   ), 0.5 , 0, 0.25, 0.95) +
        draw_plot(myBoxplotFunc(modelTable, "coverage"    , "darkgreen"), 0.75, 0, 0.25, 0.95) +
        draw_label("All genes, regression analysis", x = 0.5, y = 0.975)
)

modelTable <- inner_join(modelTable, refTable, by = "gene_id")

geneTypes <- table(modelTable$gene_type) %>% sort(decreasing = TRUE) %>% subset(., . >= 1) %>% names

t0 <- Sys.time()
pdf("../../../plots/slopes_tss_allGeneType.pdf", height = 4, width = 14)
p <- ggdraw() +
  draw_plot(myBoxplotFunc(modelTable, "CpG_sites"   , "red"      ), 0,    0, 0.25, 0.95) +
  draw_plot(myBoxplotFunc(modelTable, "mCpG_ratio"  , "blue"     ), 0.25, 0, 0.25, 0.95) +
  draw_plot(myBoxplotFunc(modelTable, "mCpG_density", "purple"   ), 0.5 , 0, 0.25, 0.95) +
  draw_plot(myBoxplotFunc(modelTable, "coverage"    , "darkgreen"), 0.75, 0, 0.25, 0.95) +
  draw_label("All genes, regression analysis", x = 0.5, y = 0.975)
print(p)
for (i in geneTypes) {
  myTable <- filter(modelTable, gene_type == i)
  p <- ggdraw() +
      draw_plot(myBoxplotFunc(myTable, "CpG_sites"   , "red"      ), 0,    0, 0.25, 0.95) +
      draw_plot(myBoxplotFunc(myTable, "mCpG_ratio"  , "blue"     ), 0.25, 0, 0.25, 0.95) +
      draw_plot(myBoxplotFunc(myTable, "mCpG_density", "purple"   ), 0.5 , 0, 0.25, 0.95) +
      draw_plot(myBoxplotFunc(myTable, "coverage"    , "darkgreen"), 0.75, 0, 0.25, 0.95) +
      draw_label(paste0(i, ", regression analysis"), x = 0.5, y = 0.975)
  print(p)
}
dev.off()
Sys.time() - t0



# p <- plot_grid(
#     ggplot(data=modelTable, aes(x = cv_exp         )) + geom_density(),
#     ggplot(data=modelTable, aes(x = cv_CpG_sites   )) + geom_density(),
#     ggplot(data=modelTable, aes(x = cv_mCpG_ratio  )) + geom_density(),
#     ggplot(data=modelTable, aes(x = cv_mCpG_density)) + geom_density(),
#     ggplot(data=modelTable, aes(x = cv_coverage    )) + geom_density()
# )
#
# p <- plot_grid(
#     ggplot(data=modelTable, aes(x = r2_CpG_sites   )) + geom_density(),
#     ggplot(data=modelTable, aes(x = r2_mCpG_ratio  )) + geom_density(),
#     ggplot(data=modelTable, aes(x = r2_mCpG_density)) + geom_density(),
#     ggplot(data=modelTable, aes(x = r2_coverage    )) + geom_density()
# )
#
# ggplot(data=modelTable, aes(x = sl_mCpG_ratio)) + geom_density(from=-500, to = 500)
#
# modelTable %>% mutate(bin = ntile(cv_exp, 10)) %>%
#     ggplot(aes(factor(bin), sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-10, 10))
#
# pall <- modelTable %>% mutate(bin = ntile(r2_mCpG_ratio, 10)) %>%
#     ggplot(aes(factor(bin), sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# pfiltered <- modelTable %>% filter(cv_exp >= 1) %>% mutate(bin = ntile(r2_mCpG_ratio, 10)) %>%
#     ggplot(aes(factor(bin), sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# pfiltered2 <- modelTable %>% filter(cv_mCpG_ratio >= 1) %>% mutate(bin = ntile(r2_mCpG_ratio, 10)) %>%
#     ggplot(aes(factor(bin), sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# pfiltered3 <- modelTable %>% filter(cv_exp >= 1 & cv_mCpG_ratio >= 1) %>% mutate(bin = ntile(r2_mCpG_ratio, 10)) %>%
#     ggplot(aes(factor(bin), sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# plot_grid(pall, pfiltered, pfiltered2, pfiltered3)




# nice, seems to work
# pall <- modelTable %>% mutate(bin = cvFilter(r2_mCpG_ratio)) %>%
#     ggplot(aes(bin, sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# pfiltered <- modelTable %>% filter(cv_exp >= 1) %>% mutate(bin = cvFilter(r2_mCpG_ratio)) %>%
#     ggplot(aes(bin, sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# pfiltered2 <- modelTable %>% filter(cv_mCpG_ratio >= 1) %>% mutate(bin = cvFilter(r2_mCpG_ratio)) %>%
#     ggplot(aes(bin, sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# pfiltered3 <- modelTable %>% filter(cv_exp >= 1 & cv_mCpG_ratio >= 1) %>% mutate(bin = cvFilter(r2_mCpG_ratio)) %>%
#     ggplot(aes(bin, sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-2, 2)) +
#     geom_hline(yintercept = 0)
#
# plot_grid(pall, pfiltered, pfiltered2, pfiltered3)
#



myBoxplotFunc(modelTable, "CpG_sites")

p_CpG_sites <-
    modelTable %>% mutate(bin = cvFilter(r2_CpG_sites)) %>%
    ggplot(aes(bin, sl_CpG_sites, group = bin)) + geom_boxplot() + coord_cartesian(ylim = c(-5, 5)) +
    geom_hline(yintercept = 0) + scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                                                    labels = letters[1:])

p_mCpG_ratio   <- modelTable %>% mutate(bin = cvFilter(r2_mCpG_ratio)) %>%
    ggplot(aes(bin, sl_mCpG_ratio)) + geom_boxplot() + coord_cartesian(ylim = c(-5, 5)) +
    geom_hline(yintercept = 0)
p_mCpG_density <- modelTable %>% mutate(bin = cvFilter(r2_mCpG_density)) %>%
    ggplot(aes(bin, sl_mCpG_density)) + geom_boxplot() + coord_cartesian(ylim = c(-5, 5)) +
    geom_hline(yintercept = 0)
p_coverage     <- modelTable %>% mutate(bin = cvFilter(r2_coverage)) %>%
    ggplot(aes(bin, sl_coverage)) + geom_boxplot() + coord_cartesian(ylim = c(-5, 5)) +
    geom_hline(yintercept = 0)

plot_grid(
    p_CpG_sites, p_mCpG_ratio, p_mCpG_density, p_coverage
)

