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
[`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/save_pq.md)
to write the phyloseq object in all 3 versions (one table for each slot,
a file merging each slot and an Rdata file).

``` r
save_pq(data_fungi)
```

### Import

To import a Rdata file, just use
[`load()`](https://rdrr.io/r/base/load.html) base function. In order to
import phyloseq object from a folder create using
[`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/write_pq.md)
or
[`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/save_pq.md),
please use the function
[`read_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/read_pq.md).

``` r
d <- read_pq(path = "fungi_phyloseq")
```

### Tracking sequences, clusters and samples

In bioinformatic pipeline, we often need to track the number of samples,
sequences and clusters across step in the pipeline. MiscMetabar propose
two utilities to achieve this goal :
[`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)
and a version to compute value per samples :
[`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow_samples.md).
The function
[`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)
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
#> R version 4.5.1 (2025-06-13)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Kali GNU/Linux Rolling
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.29.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=fr_FR.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=fr_FR.UTF-8        LC_COLLATE=fr_FR.UTF-8    
#>  [5] LC_MONETARY=fr_FR.UTF-8    LC_MESSAGES=fr_FR.UTF-8   
#>  [7] LC_PAPER=fr_FR.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=fr_FR.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Paris
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] MiscMetabar_0.14.4 purrr_1.1.0        dplyr_1.1.4        dada2_1.36.0      
#> [5] Rcpp_1.1.0         ggplot2_4.0.0      phyloseq_1.52.0   
#> 
#> loaded via a namespace (and not attached):
#>   [1] bitops_1.0-9                pbapply_1.7-4              
#>   [3] deldir_2.0-4                permute_0.9-8              
#>   [5] rlang_1.1.6                 magrittr_2.0.4             
#>   [7] ade4_1.7-23                 matrixStats_1.5.0          
#>   [9] compiler_4.5.1              mgcv_1.9-3                 
#>  [11] png_0.1-8                   systemfonts_1.2.3          
#>  [13] vctrs_0.6.5                 reshape2_1.4.4             
#>  [15] stringr_1.5.2               pwalign_1.4.0              
#>  [17] pkgconfig_2.0.3             crayon_1.5.3               
#>  [19] fastmap_1.2.0               XVector_0.48.0             
#>  [21] Rsamtools_2.24.1            rmarkdown_2.29             
#>  [23] UCSC.utils_1.4.0            ragg_1.5.0                 
#>  [25] xfun_0.53                   cachem_1.1.0               
#>  [27] GenomeInfoDb_1.44.3         jsonlite_2.0.0             
#>  [29] biomformat_1.36.0           rhdf5filters_1.20.0        
#>  [31] DelayedArray_0.34.1         Rhdf5lib_1.30.0            
#>  [33] BiocParallel_1.42.2         jpeg_0.1-11                
#>  [35] parallel_4.5.1              cluster_2.1.8.1            
#>  [37] R6_2.6.1                    bslib_0.9.0                
#>  [39] stringi_1.8.7               RColorBrewer_1.1-3         
#>  [41] GenomicRanges_1.60.0        jquerylib_0.1.4            
#>  [43] SummarizedExperiment_1.38.1 iterators_1.0.14           
#>  [45] knitr_1.50                  IRanges_2.42.0             
#>  [47] Matrix_1.7-4                splines_4.5.1              
#>  [49] igraph_2.1.4                tidyselect_1.2.1           
#>  [51] abind_1.4-8                 yaml_2.3.10                
#>  [53] vegan_2.7-1                 codetools_0.2-20           
#>  [55] hwriter_1.3.2.1             lattice_0.22-7             
#>  [57] tibble_3.3.0                plyr_1.8.9                 
#>  [59] Biobase_2.68.0              withr_3.0.2                
#>  [61] ShortRead_1.66.0            S7_0.2.0                   
#>  [63] evaluate_1.0.5              desc_1.4.3                 
#>  [65] survival_3.8-3              RcppParallel_5.1.11-1      
#>  [67] Biostrings_2.76.0           pillar_1.11.1              
#>  [69] MatrixGenerics_1.20.0       foreach_1.5.2              
#>  [71] stats4_4.5.1                generics_0.1.4             
#>  [73] S4Vectors_0.46.0            scales_1.4.0               
#>  [75] glue_1.8.0                  tools_4.5.1                
#>  [77] interp_1.1-6                data.table_1.17.8          
#>  [79] GenomicAlignments_1.44.0    fs_1.6.6                   
#>  [81] rhdf5_2.52.1                grid_4.5.1                 
#>  [83] ape_5.8-1                   latticeExtra_0.6-31        
#>  [85] nlme_3.1-168                GenomeInfoDbData_1.2.14    
#>  [87] cli_3.6.5                   textshaping_1.0.3          
#>  [89] S4Arrays_1.8.1              gtable_0.3.6               
#>  [91] sass_0.4.10                 digest_0.6.37              
#>  [93] BiocGenerics_0.54.0         SparseArray_1.8.1          
#>  [95] htmlwidgets_1.6.4           farver_2.1.2               
#>  [97] htmltools_0.5.8.1           pkgdown_2.1.3              
#>  [99] multtest_2.64.0             lifecycle_1.0.4            
#> [101] httr_1.4.7                  MASS_7.3-65
```
