library(here)
library(tidyverse)
library(future.apply); plan(multiprocess(workers = 12))
library(rtracklayer)

# loading data ----------------
load(here("data", "rnaseq", "salmon_quant_transcripts_genes.RData"))

gencode <- import(here("data", "annotation", "gencode.v24.annotation.hg19.gtf"), format = "gtf") %>%
    as_tibble()

gene       <- filter(gencode, type == "gene"       & gene_type == "protein_coding")
transcript <- filter(gencode, type == "transcript" & gene_type == "protein_coding")
exon       <- filter(gencode, type == "exon"       & gene_type == "protein_coding")

getMiddleExonsFor <- function(myGene, transcript, exon) {
    myTranscripts <- filter(transcript, gene_id == myGene) %>% arrange(width)
    myExons <- filter(exon, transcript_id == myTranscripts$transcript_id[1])
    # trust gencode sorting, may need some additional sorting?
    return(myExons[c(-1,-nrow(myExons)), ])
}

t0 <- Sys.time()
myExons <- bind_rows(
    parallel::mclapply(
        gene$gene_id,
        function(x) getMiddleExonsFor(x, transcript, exon),
        mc.cores = 10
    )
)
Sys.time() - t0 # 1.2 minutes

my_exons <- dplyr::filter(myExons, !(seqnames %in% c("chrX", "chrY"))) %>%
    select(seqnames, start, end, strand) %>%
    mutate(name = paste0(seqnames, ":", start, "-", end, ";", strand)) %>%
    distinct()

fe <- new.env()
fe$make_relationship_table <- function(exonline, exon_bed, annot) {
    genes <- filter(
        annot,
        seqnames == exon_bed$seqnames[exonline] & start == exon_bed$start[exonline] & end == exon_bed$end[exonline] & strand ==  exon_bed$strand[exonline]
    )$gene_id %>%
        unique()

    transcript <- map_chr(
        genes,
        function(x) {
            filter(
                annot,
                seqnames == exon_bed$seqnames[exonline] & start == exon_bed$start[exonline] & end == exon_bed$end[exonline] & strand ==  exon_bed$strand[exonline] & gene_id == x
            )$transcript_id %>%
                unique() %>%
                sub(".[[:digit:]]+$", "", .) %>%
                paste(., collapse = ", ")
        }
    )

    if(length(transcript) != length(genes)) stop(paste("Gene / transcripts length discrepency for line:", exonline))

    cbind(
        exon_bed[exonline, ],
        transcripts = transcript,
        gene = genes
    )
}

fe$make_relationship_table(52, my_exons, exon)

t0 <- Sys.time() # 12 minutes
rel_tab_exon <- parallel::mclapply(
    seq_len(nrow(my_exons)),
    fe$make_relationship_table,
    my_exons, exon,
    mc.cores = 10
) %>% bind_rows() %>% as_tibble()
Sys.time() - t0

rel_tab_exon <- distinct(rel_tab_exon) %>%
    mutate(gene = sub(".[[:digit:]]+$", "", gene))

write_tsv(rel_tab_exon, path = here("data", "annotation", "relationship_exon_transcripts_genes_gencodev22_hg19.tsv"))
