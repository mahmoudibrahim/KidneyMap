Create Single-Cell Experiment for PDGFRb+ dataset
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Setting the environment

### Internal variables

``` r
set.seed(1234)
OUTDIR <- "./mouse_PDGFRBpositive_01_create_sce_output/";
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
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    ##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    ##     setdiff, sort, table, tapply, union, unique, unsplit, which,
    ##     which.max, which.min

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

    ## Warning: package 'matrixStats' was built under R version 3.6.2

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: BiocParallel

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
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

``` r
library("scater")
```

    ## Loading required package: ggplot2

## Load data

``` r
dat <- Matrix::readMM("../data/SmartSeq/log_expression.mtx")
rowDat <- read.table("../data/SmartSeq/log_expression_rowData.txt",sep="\t",header=TRUE)
colDat <- read.table("../data/SmartSeq/log_expression_colData.txt", sep="\t",header=TRUE)

# Genes
rownames(dat) <- rowDat$x
rownames(rowDat) <- rowDat$x

# Cells
colnames(dat) <- paste0("cell",1:ncol(dat))
rownames(colDat) <- paste0("cell",1:ncol(dat))

rm(rowDat)
## In the end we do not use single-cell experiments but multi-assay
# # Create single-cell experiment object
# dat <- SingleCellExperiment(assays = list(logcounts = dat),
#               colData=colDat,
#               rowData=rowDat)
# rm(dat, rowDat, colDat)
```

Summary of cell metadata

``` r
summary(as.data.frame(colDat))
```

    ##    Annotation.Level.1                            Annotation.Level.2
    ##  Epithelial : 68      (Myo)fibroblast                     : 76     
    ##  Mesenchymal:884      Injured Vascular Smooth Muscle Cells:265     
    ##                       Mesangial Cells                     : 74     
    ##                       Parietal Epithelial Cells           : 68     
    ##                       Pericytes                           :113     
    ##                       Vascular Smooth Muscle Cells        :356     
    ##                                                                    
    ##                               Annotation.Level.3      Time.point 
    ##  Vascular Smooth Muscle Cells 1        :172      Day 10 UUO:315  
    ##  Pericytes                             :113      Day 2 UUO :320  
    ##  Injured Vascular Smooth Muscle Cells 1:112      Uninjured :317  
    ##  Renin Producing Smooth Muscle Cells   :101                      
    ##  Vascular Smooth Muscle Cells 2        : 83                      
    ##  Injured Vascular Smooth Muscle Cells 2: 77                      
    ##  (Other)                               :294                      
    ##  Core.Matrisome.Expression.Score
    ##  Min.   :0.4518                 
    ##  1st Qu.:1.1200                 
    ##  Median :1.3231                 
    ##  Mean   :1.5617                 
    ##  3rd Qu.:1.6385                 
    ##  Max.   :5.9710                 
    ## 

## Functional characterization

### Inference of Transcription Factor activities

``` r
# load regulons
df2regulon <- function(df, regulator_name="tf") {
  regulon = df %>% split(.[regulator_name]) %>% map(function(dat) {
    targets = setNames(dat$mor, dat$target)
    likelihood = dat$likelihood
    list(tfmode = targets, likelihood = likelihood)
  })
  return(regulon)
}

regulon.df <- read.table("../data/Prior/dorothea_regulon_mouse_v1.csv", sep=",", header=TRUE, stringsAsFactors = FALSE)
regulon.df <- regulon.df[regulon.df$confidence %in% c("A","B","C"),]
regul <- df2regulon(df=regulon.df)

# Calculate TF activities
TF <- viper(eset = as.matrix(dat), regulon = regul,
              nes = T, method = "scale", minsize = 4,
              eset.filter = F, adaptive.size = F, verbose=FALSE)
```

### Inference of Pathway activities

``` r
progeny.mat <- read.table("../data/Prior/progeny_matrix_mouse_v1.txt",sep=",",header=TRUE)
rownames(progeny.mat) <- progeny.mat$X
progeny.mat <- progeny.mat[which(colnames(progeny.mat)!="X")]
progeny.mat <- as.matrix(progeny.mat)

common <- intersect(rownames(dat), rownames(progeny.mat))
  
prog <- t(as.matrix(dat[common,])) %*% progeny.mat[common,]
rn <- rownames(prog)
prog <- apply(prog,2,scale)
rownames(prog) <- rn
prog <- t(prog)
  
stopifnot(colnames(dat) == colnames(prog))
```

## Create a Single-Cell Experiment

``` r
colDat$Time.point <- factor(as.character(colDat$Time.point),
             levels=c("Uninjured","Day 2 UUO", "Day 10 UUO"))
sce <- SingleCellExperiment(assays=list("logcounts"=dat),
                colData=colDat,
                altExps=list("pathway"=SummarizedExperiment(assays=list("pathway"=prog)),
                     "TF"=SummarizedExperiment(assays=list("TF"=TF)))
                )

## Alternative: MultiAssayExperiment
# library("MultiAssayExperiment")
# mdat <- MultiAssayExperiment(experiments=list("logcounts"=dat,
#                         "pathways"=prog,
#                         "TF"=TF),
#                colData=colDat)
# rm(dat,prog,TF,colDat)
# 
# #NOTE: How to access?
```

## Save data

``` r
saveRDS(sce, file=paste0(OUTDIR,"/sce.rds"))
```

## SessionInfo

``` r
sessionInfo()
```

    ## R version 3.6.1 (2019-07-05)
    ## Platform: x86_64-conda_cos6-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /home/jperales/miniconda3/envs/kidfib/lib/R/lib/libRblas.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] scater_1.14.0               ggplot2_3.1.1              
    ##  [3] scran_1.14.1                SingleCellExperiment_1.8.0 
    ##  [5] SummarizedExperiment_1.16.0 DelayedArray_0.12.0        
    ##  [7] BiocParallel_1.20.0         matrixStats_0.56.0         
    ##  [9] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
    ## [11] IRanges_2.20.0              S4Vectors_0.24.0           
    ## [13] purrr_0.3.2                 dplyr_0.8.0.1              
    ## [15] viper_1.20.0                Biobase_2.46.0             
    ## [17] BiocGenerics_0.32.0         Matrix_1.2-17              
    ## [19] rmarkdown_1.12              nvimcom_0.9-82             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] viridis_0.5.1            mixtools_1.2.0          
    ##  [3] edgeR_3.28.0             BiocSingular_1.2.0      
    ##  [5] viridisLite_0.3.0        splines_3.6.1           
    ##  [7] DelayedMatrixStats_1.8.0 assertthat_0.2.1        
    ##  [9] statmod_1.4.34           dqrng_0.2.1             
    ## [11] GenomeInfoDbData_1.2.2   vipor_0.4.5             
    ## [13] yaml_2.2.0               pillar_1.3.1            
    ## [15] lattice_0.20-38          glue_1.3.1              
    ## [17] limma_3.42.0             digest_0.6.18           
    ## [19] XVector_0.26.0           colorspace_1.4-1        
    ## [21] htmltools_0.3.6          plyr_1.8.4              
    ## [23] pkgconfig_2.0.2          zlibbioc_1.32.0         
    ## [25] scales_1.0.0             tibble_2.1.1            
    ## [27] withr_2.1.2              lazyeval_0.2.2          
    ## [29] survival_2.44-1.1        magrittr_1.5            
    ## [31] crayon_1.3.4             evaluate_0.13           
    ## [33] MASS_7.3-51.3            segmented_1.0-0         
    ## [35] class_7.3-15             beeswarm_0.2.3          
    ## [37] tools_3.6.1              stringr_1.4.0           
    ## [39] kernlab_0.9-27           munsell_0.5.0           
    ## [41] locfit_1.5-9.4           irlba_2.3.3             
    ## [43] compiler_3.6.1           e1071_1.7-1             
    ## [45] rsvd_1.0.2               rlang_0.3.4             
    ## [47] grid_3.6.1               RCurl_1.95-4.12         
    ## [49] BiocNeighbors_1.4.0      igraph_1.2.4.1          
    ## [51] bitops_1.0-6             gtable_0.3.0            
    ## [53] R6_2.4.0                 gridExtra_2.3           
    ## [55] knitr_1.22               KernSmooth_2.23-15      
    ## [57] stringi_1.4.3            ggbeeswarm_0.6.0        
    ## [59] Rcpp_1.0.1               tidyselect_0.2.5        
    ## [61] xfun_0.6

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
