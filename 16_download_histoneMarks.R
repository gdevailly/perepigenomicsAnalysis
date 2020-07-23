setwd("/groups2/joshi_grp/guillaume/cascade/data/")

library(readr)
library(dplyr)
library(purrr)

# table made from http://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated_nosubsampling/
# 2017-05-12
histone_mod_file <- read_tsv("histone.mod.tagAligned.no.downsampling.txt") %>%
    mutate(
        weight = sub("M$", "", weight) %>% as.numeric
    )

# metadata building ----------
setwd("wgbs/roadmap/")
metadata <- read_tsv("EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")
myProms <- "../../../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.light.autosomes.bed"
roadmapExp <- list(
    pc = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.pc"),
    nc = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.nc"),
    rb = read_tsv("../../../../otherProject/ChIP_heatmap/data/roadmap/57epigenomes.RPKM.rb")
)
roadmapExp <- do.call(rbind, roadmapExp)
refTable <- read_tsv(
    "../../../../annotationData/gencode.v24.annotation.hg19.middleTSStranscript.bed",
    col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)
metadata$id[!metadata$id %in% colnames(roadmapExp)] # missing E008, E017, E021, E022, check newer data?
metadata <- filter(metadata, id %in% colnames(roadmapExp))
setwd("../../")

# weight calculation -------------
histone_mod_file <- dplyr::filter(
    histone_mod_file,
    grepl(
        paste(metadata$id, collapse = "|"),
        fileName
        )
    )

filter(histone_mod_file, grepl(".gz$", fileName))$weight %>%
    sum
# ~ 95 Go


commands <- paste0(
    "wget http://egg2.wustl.edu/roadmap/data/byFileType/alignments/consolidated_nosubsampling/",
    histone_mod_file$fileName
)
setwd("histone")

t0 <- Sys.time()
map_lgl(
    commands,
    function(x) {
        system(x)
        return(TRUE)
    }
)
Sys.time() - t0 # 24 heures (O_O)'

