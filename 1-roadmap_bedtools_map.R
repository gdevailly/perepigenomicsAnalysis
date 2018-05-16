# 2016-05-02
# download all from http://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/
mv egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/
../../../../programmes/bigWigToBedGraph E003_WGBS_FractionalMethylation.bigwig E003_WGBS_FractionalMethylation.bedgraph

bedtools makewindows -g hg19.autosoms.size -w 250 -s 100 > hg19_binned_w250_s100.bed # no X, Y, M and haplotypes
sort -i hg19_binned_w250_s100.bed > hg19_binned_w250_s100_s.bed

bedtools map -a hg19_binned_w250_s100_s.bed -b E003_WGBS_FractionalMethylation.bedgraph -c 4,4,4 -o mean,count,sum -null NA > E003_WGBS_mean_count_sum.bed




# bioconductor equivalent of bedtool map?
setwd("/mnt/ris-fas1a/linux_groups2/joshi_grp/guillaume/cascade/data/wgbs/roadmap")
library(rtracklayer)
library(dplyr)
library(readr)
library(ggplot2)

metadata <- read_tsv("EG.mnemonics.name.txt")
myFiles <- list.files("egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/")
myFiles <- myFiles[grep("*.bigwig$", myFiles)]

t0 <- Sys.time()
E003 <- import.bw(paste0("egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig/", myFiles[1]))
Sys.time() - t0 # 35 sec, 51.328.201 lines

ggplot(data.frame(score = E003$score), aes(score)) +  geom_density()
summary(E003$score)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.7800  0.8600  0.7862  0.9200  1.0000
hist(E003$score)

length(which(E003$score))


