Create Single-Cell Experiment for human PDGFRb positive dataset
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Setting the environment

### Internal variables

``` r
set.seed(1234)
OUTDIR <- "./human_PDGFRbpositive_01_create_sce_output/";
if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE);
```

### Load libraries

``` r
library("Matrix")
library("viper")
```

    ## Loading required package: Biobase

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

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

``` r
# library("gridExtra")
library("dplyr")
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library("purrr")
library("GenomeInfoDbData")
library("SingleCellExperiment")
```

    ## Loading required package: SummarizedExperiment

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:Matrix':
    ## 
    ##     expand

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     reduce

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomeInfoDb

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
library("SummarizedExperiment")
library("scran")
library("scater")
```

    ## Loading required package: ggplot2

``` r
library("progeny")
library("dorothea")
library("genesorteR")
```

## Load data

``` r
dat <- Matrix::readMM("../data/pdgfrb/pdgfrbMap_UMI_counts.mtx")
rowDat <- read.table("../data/pdgfrb/pdgfrbMap_UMI_counts_rowData.txt",sep="\t",header=TRUE, stringsAsFactors = FALSE)
colDat <- read.table("../data/pdgfrb/pdgfrbMap_UMI_counts_colData.txt", sep="\t",header=TRUE, stringsAsFactors = FALSE)
```

``` r
# Genes
rownames(dat) <- rowDat$ENSEMBL.ID
rownames(rowDat) <- rowDat$ENSEMBL.ID

# Cells
colnames(dat) <- paste0("cell",1:ncol(dat))
rownames(colDat) <- paste0("cell",1:ncol(dat))

#rm(rowDat)
```

Summary of cell metadata

``` r
summary(as.data.frame(colDat))
```

    ##  Annotation.Level.1 Annotation.Level.2 Annotation.Level.3 Kidney.Function   
    ##  Length:30270       Length:30270       Length:30270       Length:30270      
    ##  Class :character   Class :character   Class :character   Class :character  
    ##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
    ##                                                                             
    ##                                                                             
    ##                                                                             
    ##  Core.Matrisome.Expression.Score  Patient.ID       
    ##  Min.   :0.0000                  Length:30270      
    ##  1st Qu.:0.5583                  Class :character  
    ##  Median :0.9029                  Mode  :character  
    ##  Mean   :0.8509                                    
    ##  3rd Qu.:1.1349                                    
    ##  Max.   :2.1439

## Create a Single-Cell Experiment

``` r
sce <- SingleCellExperiment(assays=list("counts"=dat),
                colData=colDat,
                rowData=rowDat)
```

## Normalize data

``` r
sce = scran::computeSumFactors(sce, 
                   sizes = seq(10, 200, 20), 
                   clusters = sce$Annotation.Level.3, 
                   positive = TRUE)
sce <- logNormCounts(sce)
```

### Inference of Pathway activities

``` r
# We use the matrix instead of the sce object because currently, the wrapper of progeny for
# sce objects uses the non-log transform data, which maximizes the variability of lowly
# expressed genes in UMI counts. Thus, we prefer to use log2 data
mat <- as.matrix(logcounts(sce))
rownames(mat) <- rowData(sce)$Gene.Symbol
prog <- progeny(mat, scale=TRUE, organism="Human", top=500, perm=1)
prog <- t(prog)
stopifnot(colnames(sce) == colnames(prog))
rm(mat)
```

Add alternative experiments to SCE object

``` r
altExp(sce, "pathway") <- SummarizedExperiment(assays=list("pathway"=prog))
```

## Gene sorter

``` r
sg_level2 = sortGenes(counts(sce), 
              colData(sce)$Annotation.Level.2,
              binarizeMethod = "naive")

sg_level3 = sortGenes(counts(sce), 
              colData(sce)$Annotation.Level.3,
              binarizeMethod = "naive")

assay(sce, "binary_lvl2", withDimnames=FALSE) <- sg_level2$binary
assay(sce, "binary_lvl3", withDimnames=FALSE) <- sg_level3$binary
```

# Save data

``` r
saveRDS(sce, file=paste0(OUTDIR,"/sce.rds"))
```

## SessionInfo

``` r
sessionInfo()
```

    ## R version 4.0.2 (2020-06-22)
    ## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
    ## Running under: CentOS release 6.9 (Final)
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
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] genesorteR_0.4.3            dorothea_1.0.0             
    ##  [3] progeny_1.10.0              scater_1.16.0              
    ##  [5] ggplot2_3.3.2               scran_1.16.0               
    ##  [7] SingleCellExperiment_1.10.1 SummarizedExperiment_1.18.1
    ##  [9] DelayedArray_0.14.0         matrixStats_0.56.0         
    ## [11] GenomicRanges_1.40.0        GenomeInfoDb_1.24.0        
    ## [13] IRanges_2.22.1              S4Vectors_0.26.0           
    ## [15] GenomeInfoDbData_1.2.3      purrr_0.3.4                
    ## [17] dplyr_1.0.0                 viper_1.22.0               
    ## [19] Biobase_2.48.0              BiocGenerics_0.34.0        
    ## [21] Matrix_1.2-18              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] viridis_0.5.1             mixtools_1.2.0           
    ##  [3] tidyr_1.1.0               edgeR_3.30.0             
    ##  [5] BiocSingular_1.4.0        viridisLite_0.3.0        
    ##  [7] splines_4.0.2             DelayedMatrixStats_1.10.0
    ##  [9] statmod_1.4.34            dqrng_0.2.1              
    ## [11] vipor_0.4.5               ggrepel_0.8.2            
    ## [13] yaml_2.2.1                pillar_1.4.4             
    ## [15] lattice_0.20-41           beachmat_2.4.0           
    ## [17] glue_1.4.1                limma_3.44.1             
    ## [19] digest_0.6.25             RColorBrewer_1.1-2       
    ## [21] XVector_0.28.0            colorspace_1.4-1         
    ## [23] htmltools_0.5.0           pkgconfig_2.0.3          
    ## [25] pheatmap_1.0.12           zlibbioc_1.34.0          
    ## [27] scales_1.1.1              BiocParallel_1.22.0      
    ## [29] tibble_3.0.1              generics_0.0.2           
    ## [31] ellipsis_0.3.1            withr_2.2.0              
    ## [33] mclust_5.4.6              survival_3.2-3           
    ## [35] magrittr_1.5              crayon_1.3.4             
    ## [37] evaluate_0.14             MASS_7.3-51.6            
    ## [39] segmented_1.2-0           class_7.3-17             
    ## [41] beeswarm_0.2.3            tools_4.0.2              
    ## [43] lifecycle_0.2.0           stringr_1.4.0            
    ## [45] bcellViper_1.24.0         kernlab_0.9-29           
    ## [47] munsell_0.5.0             locfit_1.5-9.4           
    ## [49] irlba_2.3.3               compiler_4.0.2           
    ## [51] e1071_1.7-3               rsvd_1.0.3               
    ## [53] rlang_0.4.6               grid_4.0.2               
    ## [55] RCurl_1.98-1.2            BiocNeighbors_1.6.0      
    ## [57] igraph_1.2.5              bitops_1.0-6             
    ## [59] rmarkdown_2.3             gtable_0.3.0             
    ## [61] R6_2.4.1                  gridExtra_2.3            
    ## [63] knitr_1.29                KernSmooth_2.23-17       
    ## [65] stringi_1.4.6             ggbeeswarm_0.6.0         
    ## [67] Rcpp_1.0.4.6              vctrs_0.3.1              
    ## [69] tidyselect_1.1.0          xfun_0.15

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
