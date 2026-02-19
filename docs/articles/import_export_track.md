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
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Pop!_OS 24.04 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
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
#> [1] MiscMetabar_0.14.6 purrr_1.2.1        dplyr_1.2.0        dada2_1.38.0      
#> [5] Rcpp_1.1.1         ggplot2_4.0.2      phyloseq_1.54.0   
#> 
#> loaded via a namespace (and not attached):
#>   [1] ade4_1.7-23                 tidyselect_1.2.1           
#>   [3] farver_2.1.2                bitops_1.0-9               
#>   [5] Biostrings_2.78.0           S7_0.2.1                   
#>   [7] fastmap_1.2.0               GenomicAlignments_1.46.0   
#>   [9] digest_0.6.39               lifecycle_1.0.5            
#>  [11] pwalign_1.6.0               cluster_2.1.8.2            
#>  [13] survival_3.8-6              magrittr_2.0.4             
#>  [15] compiler_4.5.2              rlang_1.1.7                
#>  [17] sass_0.4.10                 tools_4.5.2                
#>  [19] igraph_2.2.1                yaml_2.3.12                
#>  [21] data.table_1.18.2.1         knitr_1.51                 
#>  [23] S4Arrays_1.10.1             htmlwidgets_1.6.4          
#>  [25] interp_1.1-6                DelayedArray_0.36.0        
#>  [27] plyr_1.8.9                  RColorBrewer_1.1-3         
#>  [29] BiocParallel_1.44.0         ShortRead_1.68.0           
#>  [31] abind_1.4-8                 hwriter_1.3.2.1            
#>  [33] withr_3.0.2                 BiocGenerics_0.56.0        
#>  [35] desc_1.4.3                  grid_4.5.2                 
#>  [37] stats4_4.5.2                latticeExtra_0.6-31        
#>  [39] multtest_2.66.0             biomformat_1.38.0          
#>  [41] Rhdf5lib_1.32.0             scales_1.4.0               
#>  [43] iterators_1.0.14            MASS_7.3-65                
#>  [45] SummarizedExperiment_1.40.0 cli_3.6.5                  
#>  [47] rmarkdown_2.30              vegan_2.7-2                
#>  [49] crayon_1.5.3                ragg_1.5.0                 
#>  [51] generics_0.1.4              otel_0.2.0                 
#>  [53] RcppParallel_5.1.11-1       reshape2_1.4.5             
#>  [55] pbapply_1.7-4               ape_5.8-1                  
#>  [57] cachem_1.1.0                rhdf5_2.54.1               
#>  [59] stringr_1.6.0               splines_4.5.2              
#>  [61] parallel_4.5.2              XVector_0.50.0             
#>  [63] matrixStats_1.5.0           vctrs_0.7.1                
#>  [65] Matrix_1.7-4                jsonlite_2.0.0             
#>  [67] IRanges_2.44.0              S4Vectors_0.48.0           
#>  [69] jpeg_0.1-11                 systemfonts_1.3.1          
#>  [71] foreach_1.5.2               jquerylib_0.1.4            
#>  [73] glue_1.8.0                  pkgdown_2.2.0              
#>  [75] codetools_0.2-20            stringi_1.8.7              
#>  [77] gtable_0.3.6                deldir_2.0-4               
#>  [79] GenomicRanges_1.62.1        tibble_3.3.1               
#>  [81] pillar_1.11.1               htmltools_0.5.9            
#>  [83] Seqinfo_1.0.0               rhdf5filters_1.22.0        
#>  [85] R6_2.6.1                    textshaping_1.0.4          
#>  [87] evaluate_1.0.5              lattice_0.22-9             
#>  [89] Biobase_2.70.0              png_0.1-8                  
#>  [91] cigarillo_1.0.0             Rsamtools_2.26.0           
#>  [93] bslib_0.10.0                SparseArray_1.10.8         
#>  [95] nlme_3.1-168                permute_0.9-10             
#>  [97] mgcv_1.9-4                  xfun_0.56                  
#>  [99] fs_1.6.6                    MatrixGenerics_1.22.0      
#> [101] pkgconfig_2.0.3
```
