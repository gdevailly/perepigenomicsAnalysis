library(here) # /home/gdevailly/work/projects/cascade
library(tidyverse)

exon_md_psis <- read_tsv(here("data", "inner_exon_psis.bed"), col_names = FALSE)
exon_md_tpm  <- read_tsv(here("data", "inner_exon_epms.bed"), col_names = FALSE)
exon_md_rela <- read_tsv(here("data/annotation/relationship_exon_transcripts_genes_gencodev29_hg19.tsv"))

load("perepigenomics/data/Rdata/geneWiseData_exonTpm_H2AK5ac.RData")
data_tpm <- byFeatureData
load("perepigenomics/data/Rdata/geneWiseData_exonPsi_H2AK5ac.RData")
data_psis <- byFeatureData
rm(byFeatureData)

gene_list <- read_tsv("perepigenomics/data/gene_list.tsv", col_types = "cccc", progress = FALSE)

# psis -----------
exon_md_psis <- mutate(exon_md_psis, name = paste0(
    X1, ":", X2, "-", X3, ";", X6
))

length(which(!(exon_md_psis$name %in% exon_md_rela$name))) # 0, cool

md_psis <- left_join(
    select(exon_md_psis, X4, name),
    exon_md_rela,
    by = "name"
) %>% rename(exon_id = X4)

length(which(!(md_psis$gene %in% gene_list$ensg))) # 15, not cool

md_psis <- left_join(
    md_psis,
    gene_list,
    by = c("gene" = "ensg")
) %>% filter(gene_type == "protein_coding")

which(duplicated(md_psis$exon_id)) %>% head()

filter(md_psis, exon_id == md_psis$exon_id[203]) %>% select(gene, symbol, gene_type)

md_psis <- select(md_psis, -gene_type)
md_psis <- filter(md_psis, exon_id %in% names(data_psis))

write_tsv(md_psis, "perepigenomics/data/exon_list_psis.tsv")

# TPM -----------
exon_md_tpm <- mutate(exon_md_tpm, name = paste0(
    X1, ":", X2, "-", X3, ";", X6
))

length(which(!(exon_md_tpm$name %in% exon_md_rela$name))) # 704

md_tpm <- left_join(
    select(exon_md_tpm, X4, name),
    exon_md_rela,
    by = "name"
) %>% rename(exon_id = X4)

length(which(!(md_tpm$gene %in% gene_list$ensg))) # 719

md_tpm <- left_join(
    md_tpm,
    gene_list,
    by = c("gene" = "ensg")
) %>% filter(gene_type == "protein_coding")

which(duplicated(md_tpm$exon_id)) %>% head()

filter(md_tpm, exon_id == md_tpm$exon_id[96]) %>% select(gene, symbol, gene_type)

md_tpm <- select(md_tpm, -gene_type)
md_tpm <- filter(md_tpm, exon_id %in% names(data_tpm))

write_tsv(md_tpm, "perepigenomics/data/exon_list_tpm.tsv")
