library(here)
library(tidyverse)

rnaseqmd <- read_tsv(here("metadata_rnaseq_read_aspra.tsv"))
full_md <- read_tsv(here("metadata_rnaseq_read.tsv"))

md <- left_join(full_md, rnaseqmd, by = c("Run" = "run_accession"))
md <- filter(md, !is.na(fastq_aspera))

check <- group_by(md, Experiment, GEO) %>%
    summarise(n_run = n_distinct(Run), n_file = n(), ratio = n_file / n_run)

table(check$ratio) # \o/
# 1  2
# 14 56

se <- filter(check, ratio == 1)$Experiment
pe <- filter(check, ratio == 2)$Experiment

length(se) + length(pe) == nrow(check) # \o/

semd <- filter(md, Experiment %in% se)
pemd <- filter(md, Experiment %in% pe)

se_commands <- map_chr(
    unique(semd$Experiment),
    function(x) {
        tmp <- filter(semd, Experiment == x)
        fastq <- strsplit(tmp$fastq_aspera, "/", fixed = TRUE) %>%
            map_chr(last) %>%
            paste0("rnaseq_reads/", .) %>%
            paste(collapse = " ")
        index <- "annotation/gencode.v29.transcripts.salmon.index"
        genes <- "annotation/gencode.v29.annotation.gff3"
        paste0(
            "salmon quant -i ",
            index,
            " -l A -r ",
            fastq,
            " -p 12 -g ",
            genes,
            " --validateMappings --seqBias --gcBias --biasSpeedSamp 5 -o salmonQuant/",
            x
        )
    }
)

write_tsv(tibble(cm = se_commands), path = file.path("~", "mnt", "genotoul_grp", "guillaume", "cascade", "salmon_se.sh"), col_names = FALSE)

pe_commands <- map_chr(
    unique(pemd$Experiment),
    function(x) {
        tmp <- filter(pemd, Experiment == x)
        fastq <- strsplit(tmp$fastq_aspera, "/", fixed = TRUE) %>%
            map_chr(last)
        fastq1 <- grep("_1.fastq.gz$", fastq, value = TRUE) %>%
            paste0("rnaseq_reads/", .) %>%
            paste(collapse = " ")
        fastq2 <- grep("_2.fastq.gz$", fastq, value = TRUE) %>%
            paste0("rnaseq_reads/", .) %>%
            paste(collapse = " ")
        if (length(fastq1) != length(fastq2)) stop("Not the same numbre of read pairs!")

        index <- "annotation/gencode.v29.transcripts.salmon.index"
        genes <- "annotation/gencode.v29.annotation.gff3"
        paste0(
            "salmon quant -i ",
            index,
            " -l A -1 ",
            fastq1,
            " -2 ",
            fastq2,
            " -p 12 -g ",
            genes,
            " --validateMappings --seqBias --gcBias --biasSpeedSamp 5 -o salmonQuant/",
            x
        )
    }
)

write_tsv(tibble(cm = pe_commands), path = file.path("~", "mnt", "genotoul_grp", "guillaume", "cascade", "salmon_pe.sh"), col_names = FALSE)


salmon quant -i annotation/gencode.v24.transcripts.salmon.index -l A \
-r rnaseq_reads/SRR020287.fastq.gz \
-p 10 -g annotation/gencode.v24.annotation.gff3 \
--seqBias --gcBias --biasSpeedSamp 5 \
-o salmonQuant/
