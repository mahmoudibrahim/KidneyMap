Get amapping table of orthologs between mouse and human
================
Javier Perales-Paton - <javier.perales@bioquant.uni-heidelberg.de>

## Setting the environment

### Internal variables

``` r
set.seed(1234)
OUTDIR <- "./00_get_ortohologs_output";
if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE);
```

### Load libraries

``` r
library(biomaRt);
```

    ## Registered S3 method overwritten by 'openssl':
    ##   method      from
    ##   print.bytes Rcpp

## Get orthologs

``` r
#--- Create the BioMart
human = useMart(host="ensembl.org",
        biomart="ENSEMBL_MART_ENSEMBL", 
        dataset="hsapiens_gene_ensembl")
mouse = useMart(host="ensembl.org", 
        biomart="ENSEMBL_MART_ENSEMBL", 
        dataset="mmusculus_gene_ensembl")

#--- Attributes of interest
attributes = c("ensembl_gene_id", 
           "mmusculus_homolog_ensembl_gene", 
           "mmusculus_homolog_perc_id_r1")
orth.mouse1 = getBM(attributes, filters="with_mmusculus_homolog",values=TRUE,
                    mart = human, uniqueRows=TRUE)
```

    ## Cache found

``` r
orth.mouse_keys <- getBM(c("hgnc_symbol","ensembl_gene_id"),
                         filter="ensembl_gene_id",values=orth.mouse1$ensembl_gene_id,
                         mart = human, uniqueRows=TRUE)
```

    ## Cache found

``` r
orth.mouse_keys2 <- getBM(c("ensembl_gene_id","mgi_symbol"),
                          filter="ensembl_gene_id",
                          values=orth.mouse1$mmusculus_homolog_ensembl_gene,
                          mart = mouse, uniqueRows=TRUE)
```

    ## Cache found

``` r
# Table of orthologs
orth.mouse <- merge(orth.mouse1, orth.mouse_keys2, 
                    by.x="mmusculus_homolog_ensembl_gene",
                    by.y="ensembl_gene_id")
orth.mouse$mgi_symbol[orth.mouse$mgi_symbol==""] <- NA
orth.mouse <- na.omit(orth.mouse)

#-  Create dictionary (list of vectors)
mmu2hsa.dic <- split(orth.mouse$ensembl_gene_id, orth.mouse$mgi_symbol)
mmu2hsa.dic <- lapply(mmu2hsa.dic, function(z) unique(z))

# Create another one for gene symbols
orth.mouse <- merge(orth.mouse1,orth.mouse_keys,by.x="ensembl_gene_id",by.y="ensembl_gene_id")
orth.mouse <- merge(orth.mouse, orth.mouse_keys2, 
                    by.x="mmusculus_homolog_ensembl_gene",
                    by.y="ensembl_gene_id")
orth.mouse$mgi_symbol[orth.mouse$mgi_symbol==""] <- NA
orth.mouse <- na.omit(orth.mouse)

#-  Create dictionary (list of vectors)
mmu2hsa <- split(orth.mouse$ensembl_gene_id, orth.mouse$mgi_symbol)
mmu2hsa <- lapply(mmu2hsa, function(z) unique(z))
```

## Save

``` r
saveRDS(mmu2hsa, file=paste0(OUTDIR,"/mmu2hsa.rds"))
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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] biomaRt_2.42.0 rmarkdown_1.12 nvimcom_0.9-82
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.1           pillar_1.3.1         compiler_3.6.1      
    ##  [4] dbplyr_1.4.0         prettyunits_1.0.2    tools_3.6.1         
    ##  [7] progress_1.2.0       digest_0.6.18        bit_1.1-14          
    ## [10] tibble_2.1.1         RSQLite_2.1.1        evaluate_0.13       
    ## [13] memoise_1.1.0        BiocFileCache_1.10.0 pkgconfig_2.0.2     
    ## [16] rlang_0.3.4          DBI_1.0.0            curl_3.3            
    ## [19] yaml_2.2.0           parallel_3.6.1       xfun_0.6            
    ## [22] dplyr_0.8.0.1        stringr_1.4.0        httr_1.4.0          
    ## [25] knitr_1.22           rappdirs_0.3.1       S4Vectors_0.24.0    
    ## [28] askpass_1.0          IRanges_2.20.0       hms_0.4.2           
    ## [31] tidyselect_0.2.5     stats4_3.6.1         bit64_0.9-7         
    ## [34] glue_1.3.1           Biobase_2.46.0       R6_2.4.0            
    ## [37] AnnotationDbi_1.48.0 XML_3.98-1.19        purrr_0.3.2         
    ## [40] blob_1.1.1           magrittr_1.5         htmltools_0.3.6     
    ## [43] BiocGenerics_0.32.0  assertthat_0.2.1     stringi_1.4.3       
    ## [46] openssl_1.3          crayon_1.3.4

``` r
{                                                                                                                                                                                                           
sink(file=paste0(OUTDIR,"/sessionInfo.txt"))
print(sessionInfo())
sink()
}
```
