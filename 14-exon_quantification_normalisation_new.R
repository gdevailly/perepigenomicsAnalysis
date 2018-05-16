library(dplyr)
library(purrr)
library(readr)
library(future); plan(multiprocess)
library(cowplot)
library(svglite)

setwd("/groups2/joshi_grp/guillaume/cascade/data/")

# failures are because presence of a left tab and the end of each line
# ... I guess. I didn't actually check each line.
geneQuant_rb %<-% read_tsv("rnaseq/roadmap/57epigenomes.RPKM.rb")
geneQuant_pc %<-% read_tsv("rnaseq/roadmap/57epigenomes.RPKM.pc")
geneQuant_nc %<-% read_tsv("rnaseq/roadmap/57epigenomes.RPKM.nc")

exonQuant_rb %<-% read_tsv("rnaseq/roadmap/57epigenomes.exn.RPKM.rb")
exonQuant_pc %<-% read_tsv("rnaseq/roadmap/57epigenomes.exon.RPKM.pc")
exonQuant_nc %<-% read_tsv("rnaseq/roadmap/57epigenomes.exon.RPKM.nc")

geneQuant %<-% bind_rows(geneQuant_rb, geneQuant_pc, geneQuant_nc)
exonQuant %<-% bind_rows(exonQuant_rb, exonQuant_pc, exonQuant_nc)

# remove duplicates, because roadmap... :-/  --------------------------
length(geneQuant$gene_id) # 52298
length(unique(geneQuant$gene_id)) # 52298
length(exonQuant$exon_location) # 307326
length(unique(exonQuant$exon_location)) # 306284

# dup <- filter(exonQuant, exon_location %in% exon_location[duplicated(exon_location)]) %>%
#     select(exon_location, everything()) %>%
#     arrange(exon_location)
# dim(dup) #  2084   59
# dim(distinct(dup)) #  1051   59

keepMaxValuesIn_exonLoc <- function(mt) {
    mt_uniq <- filter(mt, !(exon_location %in% exon_location[duplicated(exon_location)]))
    mt_dup <-  filter(mt,   exon_location %in% exon_location[duplicated(exon_location)] )
    mt_dup <- group_by(mt_dup, exon_location) %>%
        summarise_all(max)
    return(bind_rows(mt_uniq, mt_dup))
}

exonQuant <- keepMaxValuesIn_exonLoc(exonQuant)
length(exonQuant$exon_location) # 306284
length(unique(exonQuant$exon_location)) # 306284

# filter cells -----------------
metadata <- read_tsv("wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")
geneQuant <- select(geneQuant, one_of(c("gene_id", metadata$id)))
exonQuant <- select(exonQuant, one_of(c( "exon_location", "gene_id", metadata$id)))

# normalize exons ----------------------
normalizeExonFor <- function(gene, gene_RPKM = geneQuant, exon_RPKM = exonQuant) {
    g_rpkm <- filter(gene_RPKM, gene_id == gene) %>%
        select(starts_with("E", ignore.case = FALSE)) %>%
        as.matrix %>%
        as.vector
    e_rpkm <- filter(exon_RPKM, gene_id == gene)
    if(nrow(e_rpkm) >= 1) {
        exon <- e_rpkm$exon_location
        e_rpkm <- select(e_rpkm, starts_with("E", ignore.case = FALSE)) %>%
            as.matrix
        n_rpkm <- apply(e_rpkm, 1, function(x) (x+1)/(g_rpkm+1)) %>% t
        return(data.frame(
            exon_location = exon,
            gene_id = gene,
            n_rpkm,
            stringsAsFactors = FALSE
        ))
    }
}

normalizeExonFor("ENSG00000055950")[,1:5]
dim(geneQuant)
unique(geneQuant$gene_id) %>% length # \o/

t0 <- Sys.time()
normQuant <- parallel::mclapply(unique(geneQuant$gene_id), normalizeExonFor, mc.preschedule = TRUE, mc.cores = 32) %>%
    bind_rows %>%
    as_data_frame
Sys.time() - t0 # 1 minute \o/

save(normQuant, file = "Rdata/normQuant.RData")
# load("Rdata/normQuant.RData")

dim(normQuant) # 306284  35
dim(exonQuant) # 306284  35

# mean variance plot ----------------------
calc_metric <- function(table, FUN) {
    dplyr::select(table, starts_with("E", ignore.case = FALSE)) %>%
        as.matrix %>%
        apply(1, FUN)
}
# test for add_metric
# a <- data_frame(x1 = letters[1:10], E1 = 1:10, E2 = rnorm(10))
# add_metric(a, mean)

geneQuant_mean %<-% calc_metric(geneQuant, mean)
geneQuant_sd   %<-% calc_metric(geneQuant, sd  )
exonQuant_mean %<-% calc_metric(exonQuant, mean)
exonQuant_sd   %<-% calc_metric(exonQuant, sd  )
normQuant_mean %<-% calc_metric(normQuant, mean)
normQuant_sd   %<-% calc_metric(normQuant, sd  )

geneQuant_cv %<-% (geneQuant_sd / geneQuant_mean)
exonQuant_cv %<-% (exonQuant_sd / exonQuant_mean)
normQuant_cv %<-% (normQuant_sd / normQuant_mean)

geneQuant <- mutate(geneQuant, mean = geneQuant_mean, sd = geneQuant_sd, cv = geneQuant_cv) %>%
    select(gene_id, mean, sd, cv, everything())
exonQuant <- mutate(exonQuant, mean = exonQuant_mean, sd = exonQuant_sd, cv = exonQuant_cv) %>%
    select(exon_location, gene_id, mean, sd, cv, everything())
normQuant <- mutate(normQuant, mean = normQuant_mean, sd = normQuant_sd, cv = normQuant_cv) %>%
    select(exon_location, gene_id, mean, sd, cv, everything())

p_gene_sd <- ggplot(geneQuant, aes(y = log10(mean + 1), x = log10(sd+1))) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "genes")
p_gene_cv <- ggplot(geneQuant, aes(y = log10(mean + 1), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "genes")
p_exon_sd <- ggplot(exonQuant, aes(y = log10(mean + 1), x = log10(sd+1))) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "exons")
p_exon_cv <- ggplot(exonQuant, aes(y = log10(mean + 1), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "exons")
p_norm_sd <- ggplot(normQuant, aes(y = log2(mean), x = log10(sd+1))) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log2(exons+1/gene+1)", title = "normalized exons")
p_norm_cv <- ggplot(normQuant, aes(y = log2(mean), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log2(exons+1/gene+1)", title = "normalized exons")

p <- plot_grid(
    p_gene_sd, p_gene_cv,
    p_exon_sd, p_exon_cv,
    p_norm_sd, p_norm_cv,
    ncol = 2
) # slow...

save_plot(
    "../plots/genes_exons_mean_variance_plot.png",
    p,
    base_width = 5,
    base_height = 7.5
)# slow...

# gene vs exons plots ---------------------------
smallGene <- select(geneQuant, gene_id, mean)
smallExon <- select(exonQuant, exon_location, gene_id, mean)
smallNorm <- select(normQuant, exon_location, gene_id, mean)
colnames(smallGene) <- c("gene_id", "gene_rpkm")
colnames(smallExon) <- c("exon_location", "gene_id", "exon_rpkm")
colnames(smallNorm) <- c("exon_location", "gene_id", "ratio")

gene_exon %<-% full_join(smallExon, smallGene, by = "gene_id")
gene_norm %<-% full_join(smallNorm, smallGene, by = "gene_id")
exon_norm %<-% full_join(smallNorm, smallExon, by = "exon_location")

p_gene_exon <- ggplot(gene_exon, aes(x = log10(gene_rpkm+1), y = log10(exon_rpkm+1)))+ geom_point(size = 0.01) +
    geom_density2d()
p_gene_norm <- ggplot(gene_norm, aes(x = log10(gene_rpkm+1), y = log2(ratio)))+ geom_point(size = 0.01) +
    geom_density2d()
p_exon_norm <- ggplot(exon_norm, aes(x = log10(exon_rpkm+1), y = log2(ratio)))+ geom_point(size = 0.01) +
    geom_density2d()

p <- plot_grid(
    p_gene_exon,
    p_gene_norm,
    p_exon_norm,
    ncol = 3
) # slow...

save_plot(
    "../plots/genes_vs_exons_mean_fpkm.png",
    p,
    base_width = 7.5,
    base_height = 2.5
) # slow...

# selecting only middle exons, see script 9-gencode_parsing_exons.R ----------------------------
myExons <- "../../annotationData/gencode.v24.annotation.hg19.middle.exons.pc.light.autosomes.bed"
myExonsTbl <- read_tsv(myExons, col_names = FALSE)
colnames(myExonsTbl) <- c("chr", "start" ,"end", "name", "score" , "strand")
myExonsTbl$exon_location <- with(myExonsTbl, paste0(
    chr,
    ":",
    start,
    "-",
    end,
    "<",
    ifelse(strand == "+", 1, -1)
))

length(myExonsTbl$exon_location) # 25599
length(unique(myExonsTbl$exon_location)) # 25103

dup <- filter(myExonsTbl, exon_location %in% exon_location[duplicated(exon_location)]) %>%
    select(exon_location, everything()) %>%
    arrange(exon_location)

dim(dup) # 959   7
dim(distinct(dup)) # 959   7
# duplicate exons are in 2 (or more) different gene_ids...
myExonsTbl <- filter(myExonsTbl, !duplicated(exon_location))


length(myExonsTbl$exon_location) # 25103
length(exonQuant$exon_location) # 306284
length(intersect(myExonsTbl$exon_location, exonQuant$exon_location)) # 21125, good enough

innerExonQuant <- inner_join(myExonsTbl, exonQuant, by = "exon_location")
save(innerExonQuant, file = "Rdata/innerExonQuant.RData")

innerExonRatio <- inner_join(myExonsTbl, normQuant, by = "exon_location")
save(innerExonRatio, file = "Rdata/innerExonRatio.RData")

# new plots with inner exons only
p_gene_sd <- ggplot(geneQuant, aes(y = log10(mean + 1), x = log10(sd+1))) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "genes")
p_gene_cv <- ggplot(geneQuant, aes(y = log10(mean + 1), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "genes")
p_exon_sd <- ggplot(innerExonQuant, aes(y = log10(mean + 1), x = log10(sd+1))) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "exons")
p_exon_cv <- ggplot(innerExonQuant, aes(y = log10(mean + 1), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "exons")
p_norm_sd <- ggplot(innerExonRatio, aes(y = log2(mean), x = log10(sd+1))) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log2(exons+1/gene+1)", title = "normalized exons")
p_norm_cv <- ggplot(innerExonRatio, aes(y = log2(mean), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log2(exons+1/gene+1)", title = "normalized exons")

p <- plot_grid(
    p_gene_sd, p_gene_cv,
    p_exon_sd, p_exon_cv,
    p_norm_sd, p_norm_cv,
    ncol = 2
) # slow...

save_plot(
    "../plots/inner_exons_genes_exons_mean_variance_plot.png",
    p,
    base_width = 5,
    base_height = 7.5
)# slow...

smallGene <- select(geneQuant, gene_id, mean)
smallExon <- select(innerExonQuant, exon_location, gene_id, mean)
smallNorm <- select(innerExonRatio, exon_location, gene_id, mean)
colnames(smallGene) <- c("gene_id", "gene_rpkm")
colnames(smallExon) <- c("exon_location", "gene_id", "exon_rpkm")
colnames(smallNorm) <- c("exon_location", "gene_id", "ratio")

gene_exon %<-% full_join(smallExon, smallGene, by = "gene_id")
gene_norm %<-% full_join(smallNorm, smallGene, by = "gene_id")
exon_norm %<-% full_join(smallNorm, smallExon, by = "exon_location")

p_gene_exon <- ggplot(gene_exon, aes(x = log10(gene_rpkm+1), y = log10(exon_rpkm+1)))+ geom_point(size = 0.01) +
    geom_density2d()
p_gene_norm <- ggplot(gene_norm, aes(x = log10(gene_rpkm+1), y = log2(ratio)))+ geom_point(size = 0.01) +
    geom_density2d()
p_exon_norm <- ggplot(exon_norm, aes(x = log10(exon_rpkm+1), y = log2(ratio)))+ geom_point(size = 0.01) +
    geom_density2d()

p <- plot_grid(
    p_gene_exon,
    p_gene_norm,
    p_exon_norm,
    ncol = 3
) # slow...

save_plot(
    "../plots/inner_exon_genes_vs_exons_mean_fpkm.png",
    p,
    base_width = 7.5,
    base_height = 2.5
) # slow...







