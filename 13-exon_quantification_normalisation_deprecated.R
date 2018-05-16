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

# dim(geneQuant) # [1]  52298    58
# dim(exonQuant) # [1] 307326    59

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
        n_rpkm <- apply(e_rpkm, 1, function(x) x/g_rpkm) %>% t
        return(data.frame(
            exon_location = exon,
            gene_id = gene,
            n_rpkm,
            stringsAsFactors = FALSE
        ))
    }
}

# normalizeExonFor("ENSG00000055950")
# dim(geneQuant)
# unique(geneQuant$gene_id) %>% length # \o/

# t0 <- Sys.time()
# normQuant <- parallel::mclapply(unique(geneQuant$gene_id), normalizeExonFor, mc.preschedule = TRUE, mc.cores = 32) %>%
#     bind_rows %>%
#     as_data_frame
# Sys.time() - t0 # 1 minute \o/

# save(normQuant, file = "Rdata/normQuant.RData")
load("Rdata/normQuant.RData")

dim(normQuant) # 307326
dim(exonQuant) # 307326

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

p_gene_sd <- ggplot(geneQuant, aes(y = log10(mean + 1), x = sd)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "genes")
p_gene_cv <- ggplot(geneQuant, aes(y = log10(mean + 1), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "genes")
p_exon_sd <- ggplot(exonQuant, aes(y = log10(mean + 1), x = sd)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "exons")
p_exon_cv <- ggplot(exonQuant, aes(y = log10(mean + 1), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log10(FPKM+1)", title = "exons")
p_norm_sd <- ggplot(normQuant, aes(y = log2(mean), x = sd)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log2(exons/gene)", title = "normalized exons") + coord_cartesian(ylim = c(0, 10))
p_norm_cv <- ggplot(normQuant, aes(y = log2(mean), x = cv)) + geom_point(size = 0.01) + geom_density2d() + labs(y = "log2(exons/gene)", title = "normalized exons") + coord_cartesian(ylim = c(0, 10))

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

# merging with gtf ----------------------------
myExons <- "../../annotationData/gencode.v24.annotation.hg19.middle.exons.pc.light.autosomes.bed"
myExonsTbl <- read_tsv(myExons, col_names = FALSE)
colnames(myExonsTbl) <- c("chr", "start" ,"end", "name", "score" , "strand")
myExonsTbl$roadmap_id <- with(myExonsTbl, paste0(
    chr,
    ":",
    start,
    "-",
    end,
    "<",
    ifelse(strand == "+", 1, -1)
))

filter(exonQuant, duplicated(exon_location))
filter(exonQuant, exon_location == "chr12:6991360-6991491<1") # WHY ROADMAP!!!

# filter(exonQuant, exon_location %in% exon_location[duplicated(exon_location)]) %>%
#     group_by(exon_location) %>%
#     summarise_if(is.numeric, funs(IQR(., na.rm = TRUE))) %>%
#     ggplot(aes(x = mean)) +
#     geom_density()
# seems good, just keep the higest value


keepMaxValuesIn_exonLoc <- function(mt) {
    mt_uniq <- filter(mt, !(exon_location %in% exon_location[duplicated(exon_location)]))
    mt_dup <-  filter(mt,   exon_location %in% exon_location[duplicated(exon_location)] )
    mt_dup <- group_by(mt_dup, exon_location) %>%
        summarise_all(max)
    return(bind_rows(mt_uniq, mt_dup))
}

exonQuantUniq <- keepMaxValuesIn_exonLoc(exonQuantUniq)



exonQuant_rb$exon_location %>% length # 1673
exonQuant_rb$exon_location %>% unique %>% length # 1673

exonQuant_pc$exon_location %>% length # 233480
exonQuant_pc$exon_location %>% unique %>% length # 233480

exonQuant_nc$exon_location %>% length # 72173
exonQuant_nc$exon_location %>% unique %>% length # 72173

rb_pc <- bind_rows(exonQuant_rb, exonQuant_pc)
rb_pc$exon_location %>% length # 235153
rb_pc$exon_location %>% unique %>% length # 235148, unsafe, sight.

rb_nc <- bind_rows(exonQuant_rb, exonQuant_nc)
rb_nc$exon_location %>% length # 73846
rb_nc$exon_location %>% unique %>% length # 73842, unsafe, sight.

nc_pc <- bind_rows(exonQuant_nc, exonQuant_pc)
nc_pc$exon_location %>% length # 305653
nc_pc$exon_location %>% unique %>% length # 304620, unsafe, sight.


# length(myExonsTbl$roadmap_id) # 25599
# length(unique(myExonsTbl$roadmap_id)) # 25103
# length(exonQuant$exon_location) # 307326
# length(unique(exonQuant$exon_location)) # 306284
# length(intersect(unique(myExonsTbl$roadmap_id), unique(exonQuant$exon_location))) # 21125, good enought I guess


# library(stringr)
# getExonLocationFrom <- function(myvect = exonQuant$exon_location[1:10]) {
#     str_split(string = myvect, pattern = "[:,<]") %>%
#         lapply(function(x){
#             c(
#                 x[1],
#                 str_split(x[2], "-"),
#                 x[3]
#             ) %>% unlist
#         }) %>%
#         bind_rows
# }







