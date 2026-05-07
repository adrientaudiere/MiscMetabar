# Import/export and track

``` r
library(MiscMetabar)
data(data_fungi)
data(data_fungi_sp_known)
```

### Export phyloseq object

You can export a phyloseq object to csv (and txt for phylogenetic tree)
files in a folder. It is possible to export each table into one file or
to merge all slot (except phytree) in one file (args `one_file = TRUE`).
Finally, if `rdata` is set to TRUE, a `rdata` file containing the
phyloseq object is also writed.

``` r
write_pq(data_fungi, path = "fungi_phyloseq")
write_pq(data_fungi, path = "fungi_phyloseq", one_file = TRUE)
write_pq(data_fungi, path = "fungi_phyloseq", rdata = TRUE)
```

Finally, you can use the function
[`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/save_pq.md)
to write the phyloseq object in all 3 versions (one table for each slot,
a file merging each slot and an Rdata file).

``` r
save_pq(data_fungi)
```

### Import

To import a Rdata file, just use
[`load()`](https://acclab.github.io/dabestr/reference/load.html) base
function. In order to import phyloseq object from a folder create using
[`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/write_pq.md)
or
[`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/save_pq.md),
please use the function
[`read_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/read_pq.md).

``` r
d <- read_pq(path = "fungi_phyloseq")
```

### Tracking sequences, clusters and samples

In bioinformatic pipeline, we often need to track the number of samples,
sequences and clusters across step in the pipeline. MiscMetabar propose
two utilities to achieve this goal :
[`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/track_wkflow.md)
and a version to compute value per samples :
[`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/track_wkflow_samples.md).
The function
[`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/track_wkflow.md)
can deal with (i) fastq and fastg.gz files, dada-class object,
derep-class object, matrix of samples x clusters (e.g. `otu_table`) and
phyloseq-class object.

``` r
track_wkflow(list(data_fungi, data_fungi_sp_known))
#>   nb_sequences nb_clusters nb_samples
#> 1      1839124        1420        185
#> 2      1106581         651        185
```

## Session information

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Pop!_OS 24.04 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Paris
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] MiscMetabar_0.16.1.9000 divent_0.5-3            purrr_1.2.2            
#> [4] dplyr_1.2.1             dada2_1.38.0            Rcpp_1.1.1             
#> [7] ggplot2_4.0.2           phyloseq_1.54.2        
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6                bitops_1.0-9               
#>   [3] pbapply_1.7-4               deldir_2.0-4               
#>   [5] permute_0.9-10              rlang_1.2.0                
#>   [7] magrittr_2.0.5              ade4_1.7-24                
#>   [9] otel_0.2.0                  matrixStats_1.5.0          
#>  [11] compiler_4.5.2              mgcv_1.9-4                 
#>  [13] png_0.1-9                   systemfonts_1.3.2          
#>  [15] vctrs_0.7.3                 reshape2_1.4.5             
#>  [17] stringr_1.6.0               pwalign_1.6.0              
#>  [19] pkgconfig_2.0.3             crayon_1.5.3               
#>  [21] fastmap_1.2.0               XVector_0.50.0             
#>  [23] Rsamtools_2.26.0            rmarkdown_2.31             
#>  [25] ragg_1.5.2                  xfun_0.57                  
#>  [27] cachem_1.1.0                cigarillo_1.0.0            
#>  [29] jsonlite_2.0.0              biomformat_1.38.3          
#>  [31] rhdf5filters_1.22.0         DelayedArray_0.36.1        
#>  [33] Rhdf5lib_1.32.0             BiocParallel_1.44.0        
#>  [35] jpeg_0.1-11                 parallel_4.5.2             
#>  [37] cluster_2.1.8.2             R6_2.6.1                   
#>  [39] bslib_0.10.0                stringi_1.8.7              
#>  [41] RColorBrewer_1.1-3          GenomicRanges_1.62.1       
#>  [43] jquerylib_0.1.4             Seqinfo_1.0.0              
#>  [45] SummarizedExperiment_1.40.0 iterators_1.0.14           
#>  [47] knitr_1.51                  IRanges_2.44.0             
#>  [49] Matrix_1.7-4                splines_4.5.2              
#>  [51] igraph_2.2.3                tidyselect_1.2.1           
#>  [53] abind_1.4-8                 yaml_2.3.12                
#>  [55] vegan_2.7-3                 codetools_0.2-20           
#>  [57] hwriter_1.3.2.1             lattice_0.22-9             
#>  [59] tibble_3.3.1                plyr_1.8.9                 
#>  [61] Biobase_2.70.0              withr_3.0.2                
#>  [63] ShortRead_1.68.0            S7_0.2.1                   
#>  [65] evaluate_1.0.5              desc_1.4.3                 
#>  [67] survival_3.8-6              RcppParallel_5.1.11-2      
#>  [69] Biostrings_2.78.0           pillar_1.11.1              
#>  [71] MatrixGenerics_1.22.0       foreach_1.5.2              
#>  [73] stats4_4.5.2                generics_0.1.4             
#>  [75] S4Vectors_0.48.1            scales_1.4.0               
#>  [77] glue_1.8.0                  tools_4.5.2                
#>  [79] interp_1.1-6                data.table_1.18.2.1        
#>  [81] GenomicAlignments_1.46.0    fs_2.0.1                   
#>  [83] rhdf5_2.54.1                grid_4.5.2                 
#>  [85] ape_5.8-1                   rbibutils_2.4.1            
#>  [87] latticeExtra_0.6-31         nlme_3.1-168               
#>  [89] cli_3.6.6                   textshaping_1.0.5          
#>  [91] S4Arrays_1.10.1             gtable_0.3.6               
#>  [93] sass_0.4.10                 digest_0.6.39              
#>  [95] BiocGenerics_0.56.0         SparseArray_1.10.10        
#>  [97] htmlwidgets_1.6.4           farver_2.1.2               
#>  [99] htmltools_0.5.9             pkgdown_2.2.0              
#> [101] multtest_2.66.0             lifecycle_1.0.5            
#> [103] MASS_7.3-65
```
