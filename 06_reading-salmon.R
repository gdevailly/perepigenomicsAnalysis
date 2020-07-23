library(here)
library(tidyverse)
library(future); plan(multiprocess(workers = 12))

# metadata ------------
rnaseqmd <- read_tsv(here("metadata_rnaseq_read_aspra.tsv"))
full_md <- read_tsv(here("metadata_rnaseq_read.tsv"))
md <- left_join(full_md, rnaseqmd, by = c("Run" = "run_accession"))
md <- filter(md, !is.na(fastq_aspera))
rm(rnaseqmd, full_md)

# loading transcripts --------
t0 <- Sys.time()
trans <- map(
    set_names(file.path(
        "~/mnt/genotoul_grp/guillaume/cascade/salmonQuant",
        unique(md$Experiment),
        "quant.sf"
    ), unique(md$Experiment)),
    ~read_tsv(.x, col_types = "ciddd", progress = FALSE)
)
Sys.time() - t0

# checks before fusion
map_int(trans, nrow) %>% table(useNA = "ifany")
map_int(trans, ncol) %>% table(useNA = "ifany")
map_lgl(trans, ~all(.x$Name == trans[[1]]$Name)) %>% table(useNA = "ifany")

tpm_trans <- map_dfc(trans, "TPM") %>% as.matrix()
rownames(tpm_trans) <- trans[[1]]$Name

colSums(tpm_trans) %>% summary()

# loading genes --------
t0 <- Sys.time()
genes <- map(
    set_names(file.path(
        "~/mnt/genotoul_grp/guillaume/cascade/salmonQuant",
        unique(md$Experiment),
        "quant.genes.sf"
    ), unique(md$Experiment)),
    ~read_tsv(.x, col_types = "cdddd", progress = FALSE)
)
Sys.time() - t0

# checks before fusion
map_int(genes, nrow) %>% table(useNA = "ifany")
map_int(genes, ncol) %>% table(useNA = "ifany")
map_lgl(genes, ~all(.x$Name == genes[[1]]$Name)) %>% table(useNA = "ifany")

tpm_genes <- map_dfc(genes, "TPM") %>% as.matrix()
rownames(tpm_genes) <- genes[[1]]$Name

colSums(tpm_genes) %>% summary()

# median of biological repplicates ----------------
merge_md <- select(md, Experiment, id, short, name) %>%
    distinct()
table(merge_md$id)

# side effects!
median_of <- function(mid, exp_mat) {
    cols <- filter(merge_md, id == mid)$Experiment
    if (length(cols) == 1) {
        return(exp_mat[, 1, drop = TRUE])
    }
    apply(exp_mat[, cols], 1, median)
}

t0 <- Sys.time() # 5 minutes
tpm_trans_med <- map_dfr(
    set_names(unique(merge_md$id)),
    ~median_of(.x, tpm_trans)
)
Sys.time() - t0
tpm_trans_med <- as.matrix(tpm_trans_med)
rownames(tpm_trans_med) <- rownames(tpm_trans)

t0 <- Sys.time() # 1 minutes 30s
tpm_genes_med <- map_dfr(
    set_names(unique(merge_md$id)),
    ~median_of(.x, tpm_genes)
)
Sys.time() - t0
tpm_genes_med <- as.matrix(tpm_genes_med)
rownames(tpm_genes_med) <- rownames(tpm_genes)

save(tpm_trans_med, tpm_genes_med, file = here("data", "rnaseq", "salmon_quant_transcripts_genes.RData"))
load(here("data", "rnaseq", "salmon_quant_transcripts_genes.RData"))
