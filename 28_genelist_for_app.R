library(here)
library(tidyverse)

genes <- read_tsv(here("data", "annotation", "gencode.v29.annotation.hg19.middleTSS.bed"), col_names = FALSE)
genes2 <- read_tsv(here("data", "annotation", "gencode.v29.annotation.hg19.middleTSStranscript.light.autosomes.bed"), col_names = FALSE)

colnames(genes) <- c("chr", "tss", "tss+1", "ensg", "score", "strand", "gene_type", "symbol")
table(genes$chr)

colnames(genes2) <- c("chr", "start", "end", "ensg", "score", "strand")
genes2$ensg <-  sub("\\.[0-9]*", "", genes2$ensg)
genes2 <- mutate(genes2, gene_length = end - start) %>%
    mutate(length_type = case_when(
        end - start <= 1000 ~ "short",
        end - start <= 3000 ~ "intermediate",
        end - start > 3000 ~ "long",
        TRUE ~ NA_character_
    ))

load(here("perepigenomics", "data", "Rdata", "geneWiseData_tss_H2A.Z.RData"))
names(byFeatureData[1:5])

gene_list <- genes %>%
    filter(chr %in% paste0("chr", 1:22)) %>%
    select(ensg, symbol, gene_type) %>%
    mutate(ensg = sub(".[[:digit:]]+$", "", .$ensg))

gene_list <- filter(gene_list, ensg %in% names(byFeatureData))

setdiff(gene_list$ensg, names(byFeatureData)) %>% length()
setdiff(names(byFeatureData), gene_list$ensg) %>% length()

table(gene_list$gene_type) %>% sort(decreasing = TRUE) %>% sort(decreasing = TRUE) %>%
{
    x <- .
    names(x[which(x > 400)])
}

my_gene_types <- c(
    "protein_coding", "processed_pseudogene", "lincRNA", "antisense",
    "snRNA", "unprocessed_pseudogene", "misc_RNA", "miRNA", "snoRNA",
    "rRNA", "transcribed_unprocessed_pseudogene", "other"
)

gene_list <- mutate(gene_list, gene_type = if_else(gene_type %in% my_gene_types, gene_type, "other"))

gene_list <- left_join(gene_list, select(genes2, ensg, length_type), by = "ensg")
table(gene_list$length_type, useNA = "ifany") # \o/

write_tsv(gene_list, here("perepigenomics", "data", "gene_list.tsv"))

