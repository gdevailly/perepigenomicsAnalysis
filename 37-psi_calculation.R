library(here)
library(tidyverse)
library(future.apply); plan(multiprocess(workers = 12))

relat <- read_tsv(here("data", "annotation", "relationship_exon_transcripts_genes_gencodev22_hg19.tsv"))
load(here("data", "rnaseq", "salmon_quant_transcripts_genes.RData"))

rownames(tpm_genes_med) <- sub(".[[:digit:]]+$", "", rownames(tpm_genes_med))
rownames(tpm_trans_med) <- sub(".[[:digit:]]+$", "", rownames(tpm_trans_med))

plot(sort(rowMeans(tpm_genes_med)), type = "l", log = "y")
# psi not defined if tpm_genes == 0
keep <- apply(tpm_genes_med, 1, function(x) all(x > 0))
tpm_genes_med <- tpm_genes_med[keep, ]
relat <- filter(relat, gene %in% rownames(tpm_genes_med))

compute_psi <- function(line, relat, trans, genes) {
    message(line)
    my_gene <- relat$gene[line]
    gene_tpms <- genes[my_gene, , drop = TRUE]
    my_transcripts <- strsplit(relat$transcripts[line], split = ", ")[[1]]
    my_transcripts <- my_transcripts[my_transcripts %in% rownames(trans)]
    if(length(my_transcripts) == 0) return(rep(NA_real_, ncol(genes)))
    tran_tpms <- colSums(trans[my_transcripts, , drop = FALSE])
    tran_tpms / gene_tpms
}

compute_psi(42, relat, trans = tpm_trans_med, genes = tpm_genes_med)

t0 <- Sys.time() # 45 sec
psis <- parallel::mclapply(
    seq_len(nrow(relat)),
    function(x) compute_psi(x, relat, trans = tpm_trans_med, genes = tpm_genes_med),
    mc.cores = 12
) %>% do.call(rbind, .)
Sys.time() - t0

rownames(psis) <- relat$name
tokeep <- apply(psis, 1, function(x) !any(is.na(x))) # 5460
psis <- psis[tokeep, ]

summary(as.vector(psis))

p <- gather(as_tibble(psis), key = "celltype", value = "psi") %>%
    ggplot(aes(x = psi)) +
    geom_density() +
    facet_wrap(~celltype) +
    theme_bw()
ggsave("plots/psi_densities.png", p)

psis[psis > 1] <- 1 # calculation errors by Salmon ?

p <- tibble(
    exon = rownames(psis),
    mean_psi = rowMeans(psis),
    sd_psi = apply(psis, 1, sd)
) %>%
    ggplot(aes(x = sd_psi, y = mean_psi)) +
    geom_point(alpha = 0.2, size = 0.5) +
    geom_density2d(color = "magenta") +
    labs(x = "s.d.", y = "mean", title = "Exon PSI, Roadmap Epigenomics dataset") +
    theme_bw()

ggsave("plots/psi_mean_var.png", p, width = 4, height = 3.5)

bed <- tibble(
    name = rownames(psis),
    chr = strsplit(name, ":") %>% map_chr(first),
    start =  strsplit(name, ":") %>% map_chr(2) %>% strsplit("-") %>% map_chr(first),
    end   =  strsplit(name, ":") %>% map_chr(2) %>% strsplit("-") %>% map_chr(2) %>% strsplit(";") %>% map_chr(first),
    strand = strsplit(name, ";") %>% map_chr(last),
    score = 2
) %>% select(chr, start, end, name, score, strand) %>%
    mutate(name = paste0("exon_id_", seq_len(nrow(.))))

write_tsv(bed, path = "data/inner_exon_psis.bed", col_names = FALSE)

psis2 <- psis
rownames(psis2) <- paste0("exon_id_", seq_len(nrow(psis2)))
saveRDS(psis2, file = here("data", "Rdata", "innerExonPsi.rds"))
