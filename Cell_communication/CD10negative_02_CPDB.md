human CD10 negative: Run CPDB to define a landscape of cell
communication
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Setup environment

### Internal variables

``` r
set.seed(1234)
# Please change it to your local installation of CellPhoneDB
CPDB_DIR <- "~/.cpdb/releases/v2.0.0/"
```

### File structure

``` r
options(stringsAsFactors = FALSE)
# Output directory
OUTDIR <- "./CD10negative_02_CPDB_output/";
if(!dir.exists(OUTDIR)) dir.create(OUTDIR);

# We define independent output folders for lvl2 and lvl3 cell type
# annotations
OUTLVL2 <- paste0(OUTDIR,"/lvl2")
if(!dir.exists(OUTLVL2)) dir.create(OUTLVL2);
OUTLVL3 <- paste0(OUTDIR,"/lvl3")
if(!dir.exists(OUTLVL3)) dir.create(OUTLVL3);

# # Figures
# FIGDIR <- paste0(OUTDIR, "/figures/")
# knitr::opts_chunk$set(fig.path=FIGDIR)
# knitr::opts_chunk$set(dev=c('png','tiff'))
# # Data
# DATADIR <- paste0(OUTDIR, "/data/")
# if(!dir.exists(DATADIR)) dir.create(DATADIR);
```

### Load libraries

``` r
library(Matrix)
library(SingleCellExperiment)
```

    ## Loading required package: SummarizedExperiment

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     which

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
library(scran)
library(scater)
```

    ## Loading required package: ggplot2

## Load data

``` r
sce <- readRDS("./CD10negative_01_create_sce_output/sce.rds")
is(sce)
```

    ## [1] "SingleCellExperiment"       "RangedSummarizedExperiment"
    ## [3] "SummarizedExperiment"       "Vector"                    
    ## [5] "Annotated"                  "vector_OR_Vector"

``` r
dim(sce)
```

    ## [1] 23446 51849

## Extract information and map genes onto human

``` r
# take raw data and normalise it for the subset of cells
count_norm <- as.matrix(logcounts(sce))

meta_data2 <- data.frame("Cell"=rownames(colData(sce)),
                        "cell_type"=colData(sce)$Annotation.Level.2)
meta_data3 <- data.frame("Cell"=rownames(colData(sce)),
                        "cell_type"=colData(sce)$Annotation.Level.3)
```

## Save CPDB input

``` r
# LVL2
if(!file.exists(paste0(OUTLVL2,"/cellphonedb_count.txt"))) {
write.table(count_norm, paste0(OUTLVL2,"/cellphonedb_count.txt"), sep="\t", quote=F, col.names=NA)
}

if(!file.exists(paste0(OUTLVL2,"/cellphonedb_meta.txt"))) {
write.table(meta_data2, paste0(OUTLVL2,"/cellphonedb_meta.txt"), sep="\t", quote=F, row.names=F)
}
if(!file.exists(paste0(OUTLVL2,"/PBS_CellPhoneDB.sh"))) {
file.copy(from = "../templates/PBS_CellPhoneDB.sh",
      to = paste0(OUTLVL2,"/PBS_CellPhoneDB.sh"))
}
```

    ## [1] TRUE

``` r
# LVL3
if(!file.exists(paste0(OUTLVL3,"/cellphonedb_count.txt"))) {
write.table(count_norm, paste0(OUTLVL3,"/cellphonedb_count.txt"), sep="\t", quote=F, col.names=NA)
}

if(!file.exists(paste0(OUTLVL3,"/cellphonedb_meta.txt"))) {
write.table(meta_data2, paste0(OUTLVL3,"/cellphonedb_meta.txt"), sep="\t", quote=F, row.names=F)
}
if(!file.exists(paste0(OUTLVL3,"/PBS_CellPhoneDB.sh"))) {
file.copy(from = "../templates/PBS_CellPhoneDB.sh",
      to = paste0(OUTLVL3,"/PBS_CellPhoneDB.sh"))
}
```

    ## [1] TRUE

## Run CPDB

``` r
# LVL2
cmd <- paste0("qsub ",paste0(OUTLVL2,"/PBS_CellPhoneDB.sh"))
if(!dir.exists(paste0(OUTLVL2,"/out"))) {
    #system(cmd) # To be shuttled in the HPC
}
```

    ## NULL

``` r
# LVL3
cmd <- paste0("qsub ",paste0(OUTLVL3,"/PBS_CellPhoneDB.sh"))
if(!dir.exists(paste0(OUTLVL3,"/out"))) {
    #system(cmd) # To be shuttled in the HPC
}
```

    ## NULL

## SessionInfo

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
    ## Running under: CentOS release 6.10 (Final)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /net/data.isilon/ag-saez/bq_jperales/KidneyMap/envs/kidneymap01/lib/libopenblasp-r0.3.10.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.utf8       LC_NUMERIC=C             
    ##  [3] LC_TIME=en_US.utf8        LC_COLLATE=en_US.utf8    
    ##  [5] LC_MONETARY=en_US.utf8    LC_MESSAGES=en_US.utf8   
    ##  [7] LC_PAPER=en_US.utf8       LC_NAME=C                
    ##  [9] LC_ADDRESS=C              LC_TELEPHONE=C           
    ## [11] LC_MEASUREMENT=en_US.utf8 LC_IDENTIFICATION=C      
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] scater_1.16.0               ggplot2_3.3.2              
    ##  [3] scran_1.16.0                SingleCellExperiment_1.10.1
    ##  [5] SummarizedExperiment_1.18.1 DelayedArray_0.14.0        
    ##  [7] matrixStats_0.56.0          Biobase_2.48.0             
    ##  [9] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
    ## [11] IRanges_2.22.1              S4Vectors_0.26.0           
    ## [13] BiocGenerics_0.34.0         Matrix_1.2-18              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.4.6              rsvd_1.0.3               
    ##  [3] locfit_1.5-9.4            lattice_0.20-41          
    ##  [5] digest_0.6.25             R6_2.4.1                 
    ##  [7] evaluate_0.14             pillar_1.4.4             
    ##  [9] zlibbioc_1.34.0           rlang_0.4.6              
    ## [11] irlba_2.3.3               rmarkdown_2.3            
    ## [13] BiocNeighbors_1.6.0       statmod_1.4.34           
    ## [15] BiocParallel_1.22.0       stringr_1.4.0            
    ## [17] igraph_1.2.5              RCurl_1.98-1.2           
    ## [19] munsell_0.5.0             vipor_0.4.5              
    ## [21] compiler_4.0.2            BiocSingular_1.4.0       
    ## [23] xfun_0.15                 pkgconfig_2.0.3          
    ## [25] ggbeeswarm_0.6.0          htmltools_0.5.0          
    ## [27] tidyselect_1.1.0          gridExtra_2.3            
    ## [29] tibble_3.0.1              GenomeInfoDbData_1.2.3   
    ## [31] edgeR_3.30.0              viridisLite_0.3.0        
    ## [33] withr_2.2.0               crayon_1.3.4             
    ## [35] dplyr_1.0.0               bitops_1.0-6             
    ## [37] grid_4.0.2                gtable_0.3.0             
    ## [39] lifecycle_0.2.0           magrittr_1.5             
    ## [41] scales_1.1.1              dqrng_0.2.1              
    ## [43] stringi_1.4.6             XVector_0.28.0           
    ## [45] viridis_0.5.1             limma_3.44.1             
    ## [47] DelayedMatrixStats_1.10.0 ellipsis_0.3.1           
    ## [49] vctrs_0.3.1               generics_0.0.2           
    ## [51] tools_4.0.2               glue_1.4.1               
    ## [53] beeswarm_0.2.3            purrr_0.3.4              
    ## [55] yaml_2.2.1                colorspace_1.4-1         
    ## [57] knitr_1.29

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
