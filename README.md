# Scripts for processing Roadmap Epigenomics data to feed the PEREpigenomics R shiny app

## PEREpigeneomics Public Instance

A public instance of PEREpigenomics is avaiable here:
[joshiapps.cbu.uib.no/perepigenomics_app](https://joshiapps.cbu.uib.no/perepigenomics_app/)

Source code for the application is available here:
[forgemia.inra.fr/guillaume.devailly/perepigenomics_app](https://forgemia.inra.fr/guillaume.devailly/perepigenomics_app/)

## Running instructions

R and bash scripts in this reppository should be run in the order indicated in the number in their file name.
Due to the complex history of this projects, many paths are hard-coded and need to be changed manually.
Sorry for this, we will do better in the future.

R packages needed for each scripts are mentionned on top. They can come from either CRAN or Bioconductor, and can be installed using:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("{packagename}")
```

where {packagename} is replaced by the package name. For example, to install Repitools:
```r
BiocManager::install("Repitools")
```

Do not hesitate to contact us, for example by creating an issue in this reppository.

## More informations:
A preprint describing our work will be available soon.

## sessionInfo():

The following software versions were used:

```r
> sessionInfo()
R version 3.6.3 (2020-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C               LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8     LC_MONETARY=fr_FR.UTF-8   
 [6] LC_MESSAGES=fr_FR.UTF-8    LC_PAPER=fr_FR.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] vroom_1.0.2                       Repitools_1.30.0                  raster_2.8-19                     sp_1.3-1                         
 [5] viridis_0.5.1                     viridisLite_0.3.0                 seqplots_1.22.2                   plotrix_3.7-6                    
 [9] future.apply_1.3.0                future_1.14.0                     forcats_0.4.0                     stringr_1.4.0                    
[13] purrr_0.3.3                       tidyr_1.1.0                       tibble_3.0.2                      ggplot2_3.1.1                    
[17] tidyverse_1.2.1                   here_0.1                          broom_0.5.2                       cowplot_1.0.0                    
[21] BSgenome.Hsapiens.UCSC.hg19_1.4.0 BSgenome_1.52.0                   Biostrings_2.52.0                 XVector_0.24.0                   
[25] readr_1.3.1                       dplyr_1.0.0                       rtracklayer_1.44.0                GenomicRanges_1.36.0             
[29] GenomeInfoDb_1.20.0               IRanges_2.18.0                    S4Vectors_0.22.0                  BiocGenerics_0.30.0              

loaded via a namespace (and not attached):
  [1] readxl_1.3.1                backports_1.1.4             spam_2.2-2                  R.devices_2.16.0            aroma.light_3.14.0         
  [6] plyr_1.8.4                  R.rsp_0.43.1                lazyeval_0.2.2              splines_3.6.3               BiocParallel_1.18.0        
 [11] listenv_0.7.0               digest_0.6.19               htmltools_0.3.6             gdata_2.18.0                Rsolnp_1.16                
 [16] magrittr_1.5                memoise_1.1.0               aroma.apd_0.6.0             cluster_2.1.0               limma_3.40.2               
 [21] annotate_1.62.0             globals_0.12.4              modelr_0.1.4                matrixStats_0.54.0          R.utils_2.9.0              
 [26] aroma.core_3.2.0            colorspace_1.4-1            blob_1.1.1                  rvest_0.3.4                 haven_2.1.0                
 [31] crayon_1.3.4                RCurl_1.95-4.12             kohonen_3.0.8               jsonlite_1.6                genefilter_1.66.0          
 [36] survival_3.1-12             glue_1.4.1                  R.huge_0.9.0                gtable_0.3.0                zlibbioc_1.30.0            
 [41] DelayedArray_0.10.0         R.cache_0.13.0              maps_3.3.0                  scales_1.0.0                vsn_3.52.0                 
 [46] DBI_1.0.0                   edgeR_3.26.5                Rcpp_1.0.3                  xtable_1.8-4                bit_1.1-14                 
 [51] preprocessCore_1.46.0       dotCall64_1.0-0             DT_0.7                      truncnorm_1.0-8             htmlwidgets_1.3            
 [56] httr_1.4.0                  gplots_3.0.1.1              RColorBrewer_1.1-2          ellipsis_0.3.1              pkgconfig_2.0.2            
 [61] XML_3.98-1.19               R.methodsS3_1.7.1           locfit_1.5-9.1              DNAcopy_1.58.0              AnnotationDbi_1.46.0       
 [66] tidyselect_1.1.0            rlang_0.4.6                 reshape2_1.4.3              later_0.8.0                 munsell_0.5.0              
 [71] cellranger_1.1.0            tools_3.6.3                 cli_1.1.0                   generics_0.0.2              RSQLite_2.1.1              
 [76] bit64_0.9-7                 caTools_1.17.1.2            nlme_3.1-147                mime_0.7                    R.oo_1.22.0                
 [81] xml2_1.2.0                  compiler_3.6.3              rstudioapi_0.10             affyio_1.54.0               stringi_1.4.5              
 [86] gsmoothr_0.1.7              fields_9.8-3                lattice_0.20-41             Matrix_1.2-18               vctrs_0.3.1                
 [91] pillar_1.4.4                lifecycle_0.2.0             BiocManager_1.30.4          bitops_1.0-6                httpuv_1.5.1               
 [96] R.filesets_2.13.0           affy_1.62.0                 R6_2.4.0                    promises_1.0.1              KernSmooth_2.23-17         
[101] gridExtra_2.3               aroma.affymetrix_3.2.0      codetools_0.2-16            gtools_3.8.1                Ringo_1.48.0               
[106] MASS_7.3-51.6               assertthat_0.2.1            SummarizedExperiment_1.14.0 rprojroot_1.3-2             withr_2.1.2                
[111] GenomicAlignments_1.20.0    Rsamtools_2.0.0             GenomeInfoDbData_1.2.1      hms_0.4.2                   grid_3.6.3                 
[116] class_7.3-17                Cairo_1.5-10                PSCBS_0.65.0                Biobase_2.44.0              shiny_1.3.2  
```
