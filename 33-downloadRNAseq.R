# table : https://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/?display=500&search=mRNA-Seq&sort=sample
# accessed 2018-03-20
# linked by hand, errors ?

library(here)
library(rvest)
library(tidyverse)

md <- read_tsv("metadata_33samples_rnaseq.tsv")

md <- mutate(md, GEO = strsplit(GEO, ", ", fixed = TRUE)) %>%
    unnest()

get_sra <- function(gsm) {

    message(paste("Trying:", gsm))

    gsm_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm)
    geo_wp <-  read_html(gsm_url)

    hrefs <- html_nodes(geo_wp, "body") %>%
        html_nodes("a")

    href <- hrefs[which(grepl("sra", hrefs, fixed = TRUE))]

    if(length(href) != 1) {
        message("No SRA link found")
        return(NA_character_)
    }

    str_extract(href, ">[[:alnum:]]+<") %>%
        sub(">", "", .) %>%
        sub("<", "", .)

}

md$sra <- map_chr(md$GEO, get_sra)
# Trying: GSM997220
# No SRA link found -> SRA link absent in GEO...
#
# Trying: GSM916093
# No SRA link found -> SRA link absent in GEO...
#
# Trying: GSM916094
# No SRA link found -> SRA link absent in GEO...
#
# can't find them, asked for help: https://twitter.com/G_Devailly/status/981110233245372416

md <- filter(md, !is.na(sra))

paste(md$sra, collapse = ", ")
# fed into: https://www.ncbi.nlm.nih.gov/Traces/study/
# got back: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRX007165%2C+SRX014392%2C+SRX142115%2C+SRX142116%2C+SRX142111%2C+SRX142112%2C+SRX142107%2C+SRX142108%2C+SRX142109%2C+SRX142110%2C+SRX056682%2C+SRX056683%2C+SRX142113%2C+SRX142114%2C+SRX259090%2C+SRX259092%2C+SRX259060%2C+SRX259061%2C+SRX259062%2C+SRX259063%2C+SRX259089%2C+SRX259091%2C+SRX259078%2C+SRX259080%2C+SRX259079%2C+SRX259082%2C+SRX1158130%2C+SRX1158129%2C+SRX135562%2C+SRX1158046%2C+SRX1158048%2C+SRX1158056%2C+SRX1158065%2C+SRX1158047%2C+SRX1158049%2C+SRX1158066%2C+SRX1158064%2C+SRX1158062%2C+SRX190124%2C+SRX263858%2C+SRX218942%2C+SRX1158052%2C+SRX1158088%2C+SRX533052%2C+SRX190128%2C+SRX263859%2C+SRX214033%2C+SRX213996%2C+SRX214027%2C+SRX214020%2C+SRX214017%2C+SRX214029%2C+SRX214031%2C+SRX214018%2C+SRX213997%2C+SRX214028%2C+SRX214030%2C+SRX214019%2C+SRX214032%2C+SRX190132%2C+SRX263862%2C+SRX263863%2C+SRX190110%2C+SRX190136%2C+SRX190118%2C+SRX263864%2C+SRX190120%2C+SRX190138%2C+SRX263865%2C+SRX190140%2C+SRX263866%2C+SRX263867%2C+SRX190142%2C+SRX190144%2C+SRX263868%2C+SRX190114%2C+SRX190146%2C+SRX263871%2C+SRX190112%2C+SRX263869%2C+SRX263870%2C+SRX190116%2C+SRX190148%2C+SRX263872%2C+SRX263873&go=go
# saved as "SraRunTable.txt", 2018-04-03

sra_md <- read_tsv("SraRunTable.txt")

full_md <- select(sra_md, Experiment, Run) %>%
    left_join(md, by = c("Experiment" = "sra"))

write_tsv(full_md, "metadata_rnaseq_read.tsv")

full_md <- read_tsv("metadata_rnaseq_read.tsv")
full_md



head(full_md$Run)
# [1] "SRR020287" "SRR020288" "SRR020289" "SRR020290" "SRR020291" "SRR031628"

# REST request for fastq aspera path and size
t0 <- Sys.time() # 5 minutes
aspera_info <- map_dfr(
    full_md$Run,
    function(x) {
        message(x)
        onerow <- read_tsv(
            paste0(
                "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=",
                x,
                "&result=read_run&fields=run_accession,fastq_aspera,fastq_md5,fastq_bytes"
            ),
            col_types = "cccc"
        )
        if (nrow(onerow) != 1) stop("Too many rows in REST outputs")
        nrow <- length(strsplit(onerow$fastq_aspera, ";", fixed = TRUE)[[1]])
        if (nrow == 1) return(onerow)
        map_dfr(
            seq_len(nrow),
            function(y) {
                tibble(
                    run_accession = onerow$run_accession,
                    fastq_aspera  = strsplit(onerow$fastq_aspera, ";", fixed = TRUE)[[1]][y],
                    fastq_md5     = strsplit(onerow$fastq_md5   , ";", fixed = TRUE)[[1]][y],
                    fastq_bytes   = strsplit(onerow$fastq_bytes , ";", fixed = TRUE)[[1]][y]
                )
            }
        )
    }
)
Sys.time() - t0

write_tsv(aspera_info, "metadata_rnaseq_read_aspra.tsv")
aspera_info <- read_tsv("metadata_rnaseq_read_aspra.tsv")

sum(is.na(aspera_info$fastq_aspera))
# 22
# (è__é)
# clinical data, need dbGap access
aspera_info <- filter(aspera_info, !is.na(fastq_aspera))
sum(as.numeric(aspera_info$fastq_bytes)) # 1508 Go...

# generating commands
aspera_bin <- "~/.aspera/connect/bin/ascp"
aspera_shh <- "~/.aspera/connect/etc/asperaweb_id_dsa.openssh"
target_dir <- "~/mnt/genotoul_grp/guillaume/cascade/rnaseq_reads/"

download_commands <- paste0(
    aspera_bin,
    " -QT -l 300m -P33001 -i ",
    aspera_shh,
    " era-fasp@",
    aspera_info$fastq_aspera,
    " ",
    target_dir
)

grep("SRR596106", download_commands)

write_tsv(
    tibble(command = download_commands),
    path = "~/mnt/genotoul_grp/guillaume/cascade/download.sh",
    col_names = FALSE
)
# début 2018-03-04, 12h30

do2 <- download_commands[grep("SRR596106", download_commands)[1]:length(download_commands)]
write_tsv(
    tibble(command = do2),
    path = "~/mnt/genotoul_grp/guillaume/cascade/download__2.sh",
    col_names = FALSE
)


" ~/.aspera/connect/bin/ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR104/009/SRR1045589/SRR1045589_1.fastq.gz .
