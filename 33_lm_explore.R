# setwd("~/work/projects/cascade/")
# setwd("~/work/project/perepigenomics_app")
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(purrr)
library(jsonlite)
library(UpSetR)

load("data/modelData/model_geneWiseData_tss_cascade.RData", verbose = TRUE)

sym2genes <- readr::read_tsv("data/gene_list.tsv")

tss <- list.files("data/modelData/", pattern = "_tss_", full.names = TRUE)
names(tss) <- stringr::str_extract(tss, "[:alnum:]+\\.") %>%
    stringr::str_replace("\\.", "")

all_tables <- lapply(tss, function(x) {
    load(x)
    modelTable
})

cvFilter <- function(x, thresholds = c(0.25, 0.5, 0.75), values = c(0.125, 0.375, 0.625, 0.875)) {
    y <- rep(values[1], length(x))
    for(i in seq_along(thresholds)) {
        y[which(x > thresholds[i])] <- values[i + 1]
    }
    return(y)
}

# R2 > 0.05 ------------------
# topgenes <- lapply(all_tables, function(x) {
#     if(grepl("dnase", colnames(x)[3])) {
#         filter(x, r2_dnase >= 0.5) %>%
#             select(gene_id, sl = sl_dnase, r2 = r2_dnase)
#     } else {
#         filter(x, r2_mCpG_ratio >= 0.5) %>%
#             select(gene_id, sl = sl_mCpG_ratio, r2 = r2_mCpG_ratio)
#     }
# })

# map_int(topgenes, nrow)
# # cascade     Dnase       H2A   H2AK5ac H2BK120ac  H2BK12ac  H2BK15ac   H2BK5ac   H3K14ac   H3K18ac   H3K23ac   H3K27ac 
# # 459      2629     16717      8169      9031     12929     17698      8344     11175     10170      7139      1905 
# # H3K27me3  H3K36me3    H3K4ac   H3K4me1   H3K4me2   H3K4me3  H3K79me1  H3K79me2    H3K9ac   H3K9me3    H4K8ac   H4K91ac 
# # 39       188      8467       482      7449       763      7741     19735      3216        97      7734     11384 
# 
# map_dbl(topgenes, function(x) {
#     t.test(x$sl)$p.value
# })
# # cascade         Dnase           H2A       H2AK5ac     H2BK120ac      H2BK12ac      H2BK15ac       H2BK5ac       H3K14ac 
# # 1.212128e-01  2.009909e-88  3.275023e-10  1.233867e-43 8.550815e-129  4.467366e-85  5.488382e-19 2.109539e-100  4.933444e-35 
# # H3K18ac       H3K23ac       H3K27ac      H3K27me3      H3K36me3        H3K4ac       H3K4me1       H3K4me2       H3K4me3 
# # 2.809736e-171  1.305272e-36 1.258081e-103  3.654662e-02  2.037709e-16 3.477257e-103  3.649604e-12  7.192870e-03  1.224769e-34 
# # H3K79me1      H3K79me2        H3K9ac       H3K9me3        H4K8ac       H4K91ac 
# # 7.659483e-39  5.647410e-13  7.657762e-07  3.375014e-05  4.822384e-32  2.127737e-20 


# top 1000 genes ------------------
topgenes <- lapply(all_tables, function(x) {
    if(grepl("dnase", colnames(x)[3])) {
        slice_max(x, order_by = r2_dnase, n = 1000) %>%
            select(gene_id, sl = sl_dnase, r2 = r2_dnase)
    } else {
        slice_max(x, order_by = r2_mCpG_ratio, n = 1000) %>%
            select(gene_id, sl = sl_mCpG_ratio, r2 = r2_mCpG_ratio)
    }
})

map_dbl(topgenes, ~min(.$r2))
# cascade     Dnase       H2A   H2AK5ac H2BK120ac  H2BK12ac  H2BK15ac   H2BK5ac   H3K14ac   H3K18ac   H3K23ac   H3K27ac  H3K27me3 
# 0.3860235 0.7236943 0.9986152 0.9034900 0.9211538 0.9702389 0.9983331 0.9175932 0.9618331 0.9289018 0.8678691 0.6623516 0.2119458 
# H3K36me3    H3K4ac   H3K4me1   H3K4me2   H3K4me3  H3K79me1  H3K79me2    H3K9ac   H3K9me3    H4K8ac   H4K91ac 
# 0.2846244 0.9115871 0.3801055 0.8852102 0.4385988 0.9210413 0.9988855 0.7053307 0.2611697 0.8995291 0.9712508
map_dbl(topgenes, function(x) {
    t.test(x$sl)$p.value
})
# cascade        Dnase          H2A      H2AK5ac    H2BK120ac     H2BK12ac     H2BK15ac      H2BK5ac      H3K14ac      H3K18ac 
# 1.821931e-01 2.275127e-59 4.109449e-04 1.776217e-06 3.091160e-30 7.011307e-13 2.439160e-03 1.561479e-25 4.624005e-13 5.503054e-27 
# H3K23ac      H3K27ac     H3K27me3     H3K36me3       H3K4ac      H3K4me1      H3K4me2      H3K4me3     H3K79me1     H3K79me2 
# 4.037314e-11 1.390245e-61 5.885144e-12 1.468241e-67 4.885527e-35 3.670818e-28 1.572691e-03 1.270595e-48 1.284627e-18 7.203401e-04 
# H3K9ac      H3K9me3       H4K8ac      H4K91ac 
# 4.854623e-04 1.153268e-38 2.399746e-11 1.431729e-04 
map_dbl(topgenes, ~median(.$sl))
# cascade       Dnase         H2A     H2AK5ac   H2BK120ac    H2BK12ac    H2BK15ac     H2BK5ac     H3K14ac     H3K18ac     H3K23ac 
# -0.52248219  0.34352067  0.07842594  0.27459459  0.26521041  0.24336084  0.11347718  0.26127992  0.17344130  0.43060115  0.26539625 
# H3K27ac    H3K27me3    H3K36me3      H3K4ac     H3K4me1     H3K4me2     H3K4me3    H3K79me1    H3K79me2      H3K9ac     H3K9me3 
# 0.25703130  0.11689297  0.48909414  0.25183714  0.24694724  0.05792664  0.37244813  0.16307284  0.15797093  0.08523551  0.22001993 
# H4K8ac     H4K91ac 
# 0.10337689  0.09904517 

# dfp <- bind_rows(imap(topgenes, ~mutate(.x, mark = .y)))
# ggplot(dfp, aes(x = mark, y = sl)) +
#     geom_boxplot(outlier.shape = NA) +
#     coord_cartesian(ylim = c(-4, 4)) +
#     theme_bw(base_size = 14) +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))

# map(topgenes, ~.$gene_id) %>%
#     fromList() %>%
#     upset(nsets = 40)


# jsonlite::fromJSON("http://www.pantherdb.org/services/oai/pantherdb/supportedgenomes")
# # name taxon_id short_name                    version                          long_name
# # 1                  human     9606      HUMAN Reference Proteome 2020_04                       Homo sapiens
# 
# jsonlite::fromJSON("http://www.pantherdb.org/services/oai/pantherdb/supportedannotdatasets")
# # ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP
# 
# jsonlite::fromJSON("http://www.pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=ENSG00000182774,ENSG00000189043,ENSG00000163319&organism=9606&annotDataSet=ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP&enrichmentTestType=FISHER&correction=FDR")

panther_post <- lapply(
    topgenes,
    function(x) {
        paste0(
            "http://www.pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=",
            paste(x$gene_id, collapse = ","),
            "&organism=9606&annotDataSet=ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP&enrichmentTestType=FISHER&correction=FDR"
        )    
    }
)

t0 <- Sys.time()
panther_res <- lapply(
    panther_post,
    function(x) {
        res <- jsonlite::read_json(x, simplifyVector = TRUE)$results$result
        message("1 done!")
        Sys.sleep(4) # be polite with the API
        res
    }
)
Sys.time() - t0 # 4 minutes

# saveRDS(panther_res, file = "../perepigenomics_revision/pantherdb_goslimbp_top1000.RDS")
panther_res2 <- map(panther_res, function(x) {
    filter(x, fdr < 0.01)
})

map_int(panther_res2, nrow)

dfp <- bind_rows(panther_res, .id = "mark") %>%
    arrange(fdr) 
    
gocat <- unique(dfp$term$id)[1:15]

p <- filter(dfp, term$id %in% gocat) %>%
    mutate(mark = stringr::str_replace(mark, "cascade", "WGBS")) %>%
    mutate(type = case_when(
        grepl("muscle", term$label) ~ "Muscle related",
        grepl("dev", term$label) ~ "Development related",
        grepl("emb", term$label) ~ "Development related",
        grepl("pattern", term$label) ~ "Development related",
        grepl("orga", term$label) ~ "Development related",
        grepl("multicell", term$label) ~ "Development related",
        grepl("diffe", term$label) ~ "Development related",
        TRUE ~ "other"
    )) %>%
    mutate(labels =  paste0(term$id, ": ", term$label)) %>%
    mutate(labels = forcats::fct_reorder(labels, as.numeric(as.factor(type)))) %>%
    ggplot(aes(x = -log10(fdr), y = labels, fill = type)) +
    geom_col() +
    facet_wrap(~mark, nrow = 5) +
    labs(y = NULL, x = "log10(FDR)") +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave("../perepigenomics_revision/GO_terms_bp.png", plot = p, width = 8, height = 11)

# gene type freq ---------------

table(sym2genes$gene_type)

sym2genes$gene_type <- forcats::fct_infreq(sym2genes$gene_type, ) %>% forcats::fct_rev()

dfp <- map_dfr(names(all_tables), function(x) {
    all_genes <- 100 * table(sym2genes$gene_type) / nrow(sym2genes)
    if(x == "cascade") {
        r2 <- filter(all_tables[[x]], r2_mCpG_ratio >= 0.5)$gene_id
        r2 <- filter(sym2genes, ensg %in% r2)
        r2 <- 100 * table(r2$gene_type) / nrow(r2)
        top <- slice_max(all_tables[[x]], order_by = r2_mCpG_ratio, n = 1000)$gene_id
        top <- filter(sym2genes, ensg %in% top)
        top <- 100 * table(top$gene_type) / nrow(top)
    } else {
        r2 <- filter(all_tables[[x]], r2_dnase >= 0.5)$gene_id
        r2 <- filter(sym2genes, ensg %in% r2)
        r2 <- 100 * table(r2$gene_type) / nrow(r2)
        top <- slice_max(all_tables[[x]], order_by = r2_dnase, n = 1000)$gene_id
        top <- filter(sym2genes, ensg %in% top)
        top <- 100 * table(top$gene_type) / nrow(top)
    }
    bind_rows(
        tibble(type = names(all_genes), freq = all_genes, categ = "All"),
        tibble(type = names(r2), freq = r2, categ = "R2 > 0.5"),
        tibble(type = names(top), freq = top, categ = "Top 1000")
    ) %>% mutate(mark = x)
})

dfp <- mutate(dfp, freq = as.numeric(freq))

p <- filter(dfp, categ != "R2 > 0.5") %>%
    mutate(mark = stringr::str_replace(mark, "cascade", "WGBS")) %>%
    ggplot(aes(x = categ, y = freq, fill = factor(type, levels = levels(sym2genes$gene_type)))) +
    geom_col() +
    facet_wrap(~mark, nrow = 5) +
    scale_fill_manual(values = set_names(
                          rev(Polychrome::createPalette(12, seed <- c("#00a3a6", "#9dc544", "#ed6e6c"))),
                          unique(dfp$type)
                      )) +
    theme_bw() +
    labs(y = "% of genes", x = NULL, fill = "Gene type")
    
ggsave("../perepigenomics_revision/gene_types.png", plot = p, width = 8, height = 11)

# CpG islands -----------
cgi <- data.table::fread("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExtUnmasked.txt.gz")
cgi <- with(cgi, GRanges(V2, IRanges(V3, V4)))

genes <- readr::read_tsv("../perepigenomics_revision/gencode.v29.annotation.hg19.middleTSStranscript.light.autosomes.bed", col_names = F)
genesgr <- with(genes, GRanges(X1, IRanges(X2, X3), gene_id = X4))
dists <- as.data.frame(distanceToNearest(genesgr, cgi))

genes <- left_join(
    mutate(genes, X = 1:nrow(genes)),
    dists,
    by= c("X" = "queryHits")
) %>% mutate(X4 = sub("\\.[0-9]*$", "", X4))

dfp <- imap_dfr(topgenes, ~tibble(
        dist = left_join(.x, genes, by = c("gene_id" = "X4"))$distance,
        mark = .y
    )
)

dfp <- bind_rows(
    dfp, tibble(dist = genes$distance, mark = "All genes")
)
dfp <- mutate(dfp, cgi = if_else(dist < 100, TRUE, FALSE))

stats <- table(dfp$mark, dfp$cgi)
stats <- map_dbl(set_names(unique(dfp$mark)), function(x) {
    as.matrix(stats[c("All genes", x), ]) %>%
        prop.test() %>%
        .$p.value
})
stats <- tibble(mark = names(stats), pval = stats, stars = case_when(
    pval < 0.001 ~ "***",
    pval < 0.01 ~ "**",
    pval < 0.05 ~ "*",
    TRUE ~ ""
))  %>%  mutate(mark = stringr::str_replace(mark, "cascade", "WGBS"))


p <- dfp %>%
    mutate(mark = stringr::str_replace(mark, "cascade", "WGBS")) %>%
    group_by(mark) %>%
    summarise(with_cgi = sum(cgi)/length(cgi)) %>%
    ungroup() %>%
    ggplot(aes(x = reorder(mark, with_cgi), y = with_cgi, fill = if_else(mark == "All genes", "A", "B"))) +
    geom_col() +
    geom_text(data = stats, mapping = aes(x = mark, label = stars, y = 0.6 )) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("A" = "black", "B" = "grey")) +
    theme_bw() +
    theme(legend.position = "none") +
    coord_flip(ylim = c(0, 1)) +
    labs(y = "% of top 1000 correlated genes\nwith a CGI near their TSS", x = NULL)

ggsave("../perepigenomics_revision/cgi_plot.png", plot = p, width = 4, height = 6)


# old scripts -------------
sym2genes <- sym2genes %>%
    select(X4, X8) %>%
    mutate(X4 = sub("\\.[0-9]*$", "", X4)) %>%
    distinct()
g2s <- sym2genes$X8
names(g2s) <- sym2genes$X4

posgenes <- g2s[dplyr::filter(modelTable, r2_mCpG_ratio >= 0.25, sl_mCpG_ratio > 0)$gene_id] %>% sub("/", "", .)
neggenes <- g2s[dplyr::filter(modelTable, r2_mCpG_ratio >= 0.25, sl_mCpG_ratio < 0)$gene_id] %>% sub("/", "", .)
nulgenes <- g2s[dplyr::filter(modelTable, r2_mCpG_ratio < 0.25)$gene_id]                     %>% sub("/", "", .)

write.table(posgenes, file = "wgbs_tss_r2_0.25_slope_positive.txt", col.names = F, row.names = F, quote = FALSE)
write.table(neggenes, file = "wgbs_tss_r2_0.25_slope_negative.txt", col.names = F, row.names = F, quote = FALSE)
write.table(nulgenes, file = "wgbs_tss_r2_0.25_null_genes.txt"    , col.names = F, row.names = F, quote = FALSE)

# CpG islands -------------
cgi <- data.table::fread("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExtUnmasked.txt.gz")
cgi <- with(cgi, GRanges(V2, IRanges(V3, V4)))

genes <- readr::read_tsv("data/annotation/gencode.v29.annotation.hg19.middleTSS.bed", col_names = F)

genesgr <- with(genes, GRanges(X1, IRanges(X2, X3), gene_id = X4))


dists <- as.data.frame(distanceToNearest(genesgr, cgi))
genes <- left_join(
    mutate(genes, X = 1:nrow(genes)),
    dists,
    by= c("X" = "queryHits")
) %>% mutate(X4 = sub("\\.[0-9]*$", "", X4))



posgenes <- dplyr::filter(modelTable, r2_mCpG_ratio >= 0.25, sl_mCpG_ratio > 0)$gene_id
neggenes <- dplyr::filter(modelTable, r2_mCpG_ratio >= 0.25, sl_mCpG_ratio < 0)$gene_id
nulgenes <- dplyr::filter(modelTable, r2_mCpG_ratio < 0.25)$gene_id

dfp <- bind_rows(
    tibble(type = "posgenes", distance = dplyr::filter(genes, X4 %in% posgenes)$distance),
    tibble(type = "neggenes", distance = dplyr::filter(genes, X4 %in% neggenes)$distance),
    tibble(type = "nulgenes", distance = dplyr::filter(genes, X4 %in% nulgenes)$distance)
)

ggplot(dfp, aes(x = type, y = distance + 1)) +
    geom_violin() +
    geom_boxplot() +
    scale_y_log10()

