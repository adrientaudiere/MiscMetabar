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
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Pop!_OS 24.04 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
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
#> [1] MiscMetabar_0.17.0.9000 dplyr_1.2.1             ggplot2_4.0.3          
#> [4] phyloseq_1.56.0        
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6          xfun_0.58             bslib_0.11.0         
#>  [4] htmlwidgets_1.6.4     Biobase_2.72.0        lattice_0.22-9       
#>  [7] Rdpack_2.6.6          vctrs_0.7.3           tools_4.6.1          
#> [10] generics_0.1.4        biomformat_1.40.0     stats4_4.6.1         
#> [13] parallel_4.6.1        tibble_3.3.1          cluster_2.1.8.2      
#> [16] pkgconfig_2.0.3       Matrix_1.7-5          data.table_1.18.4    
#> [19] RColorBrewer_1.1-3    S7_0.2.2              desc_1.4.3           
#> [22] S4Vectors_0.50.1      RcppParallel_5.1.11-2 lifecycle_1.0.5      
#> [25] compiler_4.6.1        farver_2.1.2          stringr_1.6.0        
#> [28] textshaping_1.0.5     Biostrings_2.80.1     Seqinfo_1.2.0        
#> [31] codetools_0.2-20      permute_0.9-10        htmltools_0.5.9      
#> [34] sass_0.4.10           yaml_2.3.12           pkgdown_2.2.0        
#> [37] pillar_1.11.1         crayon_1.5.3          jquerylib_0.1.4      
#> [40] MASS_7.3-65           cachem_1.1.0          vegan_2.7-5          
#> [43] iterators_1.0.14      foreach_1.5.2         nlme_3.1-169         
#> [46] tidyselect_1.2.1      digest_0.6.39         stringi_1.8.7        
#> [49] reshape2_1.4.5        splines_4.6.1         ade4_1.7-24          
#> [52] fastmap_1.2.0         grid_4.6.1            cli_3.6.6            
#> [55] magrittr_2.0.5        survival_3.8-6        ape_5.8-1            
#> [58] withr_3.0.3           scales_1.4.0          rmarkdown_2.31       
#> [61] XVector_0.52.0        igraph_2.3.3          multtest_2.68.0      
#> [64] otel_0.2.0            ragg_1.5.2            pbapply_1.7-4        
#> [67] evaluate_1.0.5        knitr_1.51            rbibutils_2.4.1      
#> [70] IRanges_2.46.0        mgcv_1.9-4            rlang_1.2.0          
#> [73] divent_0.5-4          Rcpp_1.1.1-1.1        glue_1.8.1           
#> [76] BiocGenerics_0.58.1   jsonlite_2.0.0        R6_2.6.1             
#> [79] plyr_1.8.9            systemfonts_1.3.2     fs_2.1.0
```
