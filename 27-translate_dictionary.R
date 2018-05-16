#############################
### SET WORKING DIRECTORY ###
#############################

library(tidyverse)
available_png <- read_tsv("data/availablePng.tsv")
my_dictionary <- read_tsv("data/dictionary.tsv")
source("utils.R")

t0 <- Sys.time() # 2 minutes
translated_dic <- with(available_png,
     data_frame(
         cell_id = cell,
         cell_type = translateToHuman(cell, dict = my_dictionary),
         assay = translateToHuman(assay, dict = my_dictionary),
         feature = translateToHuman(feature, dict = my_dictionary),
         gene_type = gene_type
     )
)
Sys.time() - t0

write_tsv(arrange(translated_dic, cell_id), "data/translatedAvailablePng.tsv")
