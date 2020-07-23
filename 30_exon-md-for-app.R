library(here) # /home/gdevailly/work/projects/cascade
library(tidyverse)

exon_md_psis <- read_tsv(here("data", "inner_exon_psis.bed"), col_names = FALSE)
exon_md_fpkm <- read_tsv(here("data", "inner_exon.bed"), col_names = FALSE)
exon_md_rela <- read_tsv(here("data/annotation/relationship_exon_transcripts_genes_gencodev22_hg19.tsv"))

load("perepigenomics/data/Rdata/geneWiseData_exonFpkm_H2AK5ac.RData")
data_fpkm <- byFeatureData
load("perepigenomics/data/Rdata/geneWiseData_exonPsi_H2AK5ac.RData")
data_psis <- byFeatureData
rm(byFeatureData)

gene_list <- read_tsv("perepigenomics/data/gene_list.tsv", col_types = "ccc", progress = FALSE)

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

length(which(!(md_psis$gene %in% gene_list$ensg))) # 2031, not cool

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

# exon FPKM ---------------------

myExons <- "../../annotationData/gencode.v24.annotation.hg19.middle.exons.pc.light.autosomes.bed"
myExonsTbl <- read_tsv(myExons, col_names = FALSE)

library(future); plan(multiprocess(workers = 14))

exonQuant_rb %<-% read_tsv("data/rnaseq/roadmap/57epigenomes.exn.RPKM.rb")
exonQuant_pc %<-% read_tsv("data/rnaseq/roadmap/57epigenomes.exon.RPKM.pc")
exonQuant_nc %<-% read_tsv("data/rnaseq/roadmap/57epigenomes.exon.RPKM.nc")
exonQuant <- bind_rows(exonQuant_rb, exonQuant_pc, exonQuant_nc)
rm(exonQuant_rb, exonQuant_pc, exonQuant_nc)

md_fpkm <- select(exonQuant, exon_location, gene_id)

length(which(!(names(data_fpkm) %in% md_fpkm$exon_location))) # 0 \o/

md_fpkm <- filter(md_fpkm, exon_location %in% names(data_fpkm)) %>%
    left_join(gene_list, by = c("gene_id" = "ensg")) %>%
    select(-gene_type)

write_tsv(md_fpkm, "perepigenomics/data/exon_list_fpkm.tsv")
