setwd("/media/gdevailly/SANS TITRE/inra/cascade")

source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/6-plotingFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/11-geneWiseFunctions.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/20-functions_for_histoneMarks.R")
source("~/mnt/inra_p/projets/cascade/perepigenomicsAnalysis/29-geneWiseFunctions_hisMods.R")

library(readr)
library(parallel)
library(purrr)

# metadata building
metadata <- read_tsv("~/mnt/genotoul_grp/guillaume/cascade/data/wgbs/roadmap/EG.mnemonics.name.txt", col_names = FALSE)
colnames(metadata) <- c("id", "short", "name")

myProms <- "~/mnt/genotoul_grp/guillaume/cascade/annotation/gencode.v29.annotation.hg19.middleTSStranscript.light.autosomes.bed"

salmonExp <- read_tsv("~/mnt/genotoul_grp/guillaume/cascade/data/rnaseq/salmon_exp_genes_expressed.tsv")

refTable <- read_tsv(
    "~/mnt/genotoul_grp/guillaume/cascade/annotation/gencode.v29.annotation.hg19.middleTSStranscript.bed",
    col_names = FALSE
)
colnames(refTable) <- c("chr", "start", "end", "gene_id", "score", "strand", "gene_type", "symbol")
refTable$gene_id <-  sub("\\.[0-9]*", "", refTable$gene_id)
refTable <- filter(refTable, chr %in% paste0("chr", 1:22))

metadata$id[!metadata$id %in% colnames(salmonExp)] # missing E008, E017, E021, E022, and more
metadata <- dplyr::filter(metadata, id %in% colnames(salmonExp))

tss_table <- read_tsv(myProms, col_names = FALSE, progress = FALSE)
colnames(tss_table) <- c("chr", "start", "end", "name", "score", "strand")

TES <- "~/mnt/genotoul_grp/guillaume/cascade/annotation/gencode.v29.annotation.hg19.middleTES.light.autosomes.bed"
tes_table <- read_tsv(TES, col_names = FALSE, progress = FALSE)
colnames(tes_table) <- c("chr", "start", "end", "name", "score", "strand")
tes_table <- mutate(
    tes_table,
    new_start = if_else(strand == "+", end -1L, start),
    new_end = if_else(strand == "+", end, start + 1L)
) %>%
    dplyr::select(-start, -end) %>%
    dplyr::rename(start = new_start, end = new_end) %>%
    dplyr::select(chr, start, end, name, score, strand)


# histone ------------
hisFiles <- data_frame(
    file = list.files("histone/") %>%
        grep(".gz$", ., value = TRUE),
    name = gsub(".tagAlign.gz", "", file, fixed = TRUE)
)
hisFiles <- bind_cols(
    hisFiles,
    strsplit(hisFiles$name, "-") %>%
        map_df(~tibble(cellCode = .x[1], ChIP = .x[2]))
)
hisFiles <- filter(hisFiles, cellCode %in% colnames(salmonExp))
his_md <- left_join(
    dplyr::filter(hisFiles, !ChIP %in% c("DNase", "DNaseControl", "Input")),
    dplyr::filter(hisFiles, ChIP == "Input") %>%
        dplyr::select(file, cellCode) %>%
        dplyr::rename(input_file = file),
    by = "cellCode"
)
myHisMods <- names(table(his_md$ChIP)[table(his_md$ChIP) > 1])

# preparing data ------------------
preparDataFor <- function(thisHisMod) { # unpure

    mdForThisHisMod <- dplyr::filter(his_md, ChIP == thisHisMod)

    dataForThisHisMod <- lapply(
        seq_len(nrow(mdForThisHisMod)),
        function(i) {
            gc()
            hisModData <- prepareHisModData(i, tssAnno = tss_table, metadataTable = mdForThisHisMod, expressionTable = salmonExp) %>%
                addGeneTypeInfo(refTable)
            message(paste(mdForThisHisMod$name[i], "done!"))
            return(hisModData)
        }
    )
    names(dataForThisHisMod) <- mdForThisHisMod$name

    geneWiseData <- mclapply(
        dataForThisHisMod[[1]]$gene_id,
        function(x) extractGeneWiseDataForHistone(x, dataForThisHisMod),
        mc.cores = 8
    )
    names(geneWiseData) <- dataForThisHisMod[[1]]$gene_id

    assign(thisHisMod, geneWiseData)
    save(list = thisHisMod, file = paste0("~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_tss_", thisHisMod, ".RData"))
    message(thisHisMod)
    return(NULL)
}

t0 <- Sys.time()
walk(myHisMods, preparDataFor)
Sys.time() - t0

# preparing data TES ------------------
preparDataFor <- function(thisHisMod) { # unpure

    mdForThisHisMod <- dplyr::filter(his_md, ChIP == thisHisMod)

    dataForThisHisMod <- lapply(
        seq_len(nrow(mdForThisHisMod)),
        function(i) {
            gc()
            hisModData <- prepareHisModData(i, tssAnno = tes_table, metadataTable = mdForThisHisMod, expressionTable = salmonExp) %>%
                addGeneTypeInfo(refTable)
            message(paste(mdForThisHisMod$name[i], "done!"))
            return(hisModData)
        }
    )
    names(dataForThisHisMod) <- mdForThisHisMod$name

    geneWiseData <- mclapply(
        dataForThisHisMod[[1]]$gene_id,
        function(x) extractGeneWiseDataForHistone(x, dataForThisHisMod),
        mc.cores = 8
    )
    names(geneWiseData) <- dataForThisHisMod[[1]]$gene_id

    assign(thisHisMod, geneWiseData)
    save(list = thisHisMod, file = paste0("~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_tes_", thisHisMod, ".RData"))
    message(thisHisMod)
    return(NULL)
}

t0 <- Sys.time()
walk(myHisMods[9:22], preparDataFor)
Sys.time() - t0


# DNAse1 ----------------------------
DNAse_md <- inner_join(
    dplyr::filter(hisFiles, ChIP  == "DNase") %>%
        dplyr::select(-name, -ChIP),
    dplyr::filter(hisFiles, ChIP  == "DNaseControl") %>%
        dplyr::select(-name, -ChIP),
    by = "cellCode"
) %>% dplyr::rename(DNAse_file = file.x, Control_file = file.y)
DNAse_md <- filter(DNAse_md, cellCode %in% colnames(salmonExp))


dataForDnase <- lapply(
    DNAse_md$cellCode,
    function(cellCode) {
        gc()
        dnaseData <- prepareDNAseData(cellCode, tssAnno = tss_table, metadataTable = DNAse_md, expressionTable = salmonExp) %>%
            addGeneTypeInfo(refTable)
        message(paste(cellCode, "done!"))
        return(dnaseData)
    }
)
names(dataForDnase) <- DNAse_md$cellCode

t0 <- Sys.time()
Dnase <- mclapply(
    dataForDnase[[1]]$gene_id,
    function(x) extractGeneWiseDataForDnase(x, dataForDnase),
    mc.cores = 8
)
Sys.time() - t0
names(Dnase) <- dataForDnase[[1]]$gene_id

save(Dnase, file = "~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_tss_Dnase.RData")

## DNAse TES -----------
t0 <- Sys.time() # 36 minutes
dataForDnase <- lapply(
    DNAse_md$cellCode,
    function(cellCode) {
        dnaseData <- prepareDNAseData(cellCode, tssAnno = tes_table, metadataTable = DNAse_md, expressionTable = salmonExp) %>%
            addGeneTypeInfo(refTable)
        message(paste(cellCode, "done!"))
        return(dnaseData)
    }
)
names(dataForDnase) <- DNAse_md$cellCode
Sys.time() - t0

t0 <- Sys.time()
Dnase <- mclapply(
    dataForDnase[[1]]$gene_id,
    function(x) extractGeneWiseDataForDnase(x, dataForDnase),
    mc.cores = 8
)
Sys.time() - t0
names(Dnase) <- dataForDnase[[1]]$gene_id

save(Dnase, file = "~/mnt/genotoul/work/projects/cascade/Rdata/geneWiseData_tes_Dnase.RData")
