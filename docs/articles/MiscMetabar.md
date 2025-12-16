# Introduction

## Introduction to MiscMetabar: an R package to facilitate visualization and reproducibility in metabarcoding analysis

### Raison d’être

- Complete R packages dada2 and phyloseq
- Useful visualizations (`biplot_pq`, `circle_pq`, `upset_pq`,
  `ggvenn_pq`)
- Facilitate the use of targets package

### Quick overview

For an introduction to metabarcoding in R, Please visite the [state of
the
field](https://adrientaudiere.github.io/MiscMetabar/articles/states_of_fields_in_R.html)
vignettes. The [import, export and
track](https://adrientaudiere.github.io/MiscMetabar/articles/import_export_track.html)
vignette explains how import and export `phyloseq` object. Its also show
how to summarize useful information (number of sequences, samples and
clusters) across bioinformatic pipelines.

If you are interested in ecological metrics, see the vignettes
describing
[alpha-diversity](https://adrientaudiere.github.io/MiscMetabar/articles/alpha-div.html)
and
[beta-diversity](https://adrientaudiere.github.io/MiscMetabar/articles/beta-div.html)
analysis. The vignette [filter taxa and
samples](https://adrientaudiere.github.io/MiscMetabar/articles/filter.html)
describes some data-filtering processes using MiscMetabar and the
[reclustering](https://adrientaudiere.github.io/MiscMetabar/articles/Reclustering.html)
tutorial introduces the different way of clustering already-clustered
OTU/ASV. The vignette
[tengeler](https://adrientaudiere.github.io/MiscMetabar/articles/tengeler.html)
explore the dataset from Tengeler et al. (2020) using some MiscMetabar
functions.

For developers, I also wrote a vignette describing som [rules of
codes](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html).

#### Summarize a physeq object

``` r
library("MiscMetabar")
library("phyloseq")
data("data_fungi")
summary_plot_pq(data_fungi)
```

![](MiscMetabar_files/figure-html/example-1.png)

#### Create an interactive table of the tax_table

``` r
data("GlobalPatterns", package = "phyloseq")
tax_datatable(subset_taxa(
  GlobalPatterns,
  rowSums(GlobalPatterns@otu_table) > 100000
))
```

#### Sankey diagram of the tax_table

``` r
gp <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
sankey_pq(gp, taxa = c(1:5))
```

#### Upset plot for visualize distribution of taxa in function of samples variables

``` r
if(packageVersion("ggplot2") < "4.0.0"){
  upset_pq(gp, "SampleType", taxa = "Class")
}
```

## References

Tengeler, A.C., Dam, S.A., Wiesmann, M. et al. Gut microbiota from
persons with attention-deficit/hyperactivity disorder affects the brain
in mice. Microbiome 8, 44 (2020).
<https://doi.org/10.1186/s40168-020-00816-x>

## Session inform

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
#>   [1] bitops_1.0-9                deldir_2.0-4               
#>   [3] permute_0.9-8               rlang_1.1.6                
#>   [5] magrittr_2.0.4              ade4_1.7-23                
#>   [7] matrixStats_1.5.0           compiler_4.5.1             
#>   [9] mgcv_1.9-3                  png_0.1-8                  
#>  [11] systemfonts_1.2.3           vctrs_0.6.5                
#>  [13] reshape2_1.4.4              stringr_1.5.2              
#>  [15] pwalign_1.4.0               pkgconfig_2.0.3            
#>  [17] crayon_1.5.3                fastmap_1.2.0              
#>  [19] XVector_0.48.0              labeling_0.4.3             
#>  [21] Rsamtools_2.24.1            rmarkdown_2.29             
#>  [23] UCSC.utils_1.4.0            ragg_1.5.0                 
#>  [25] xfun_0.53                   cachem_1.1.0               
#>  [27] GenomeInfoDb_1.44.3         jsonlite_2.0.0             
#>  [29] biomformat_1.36.0           rhdf5filters_1.20.0        
#>  [31] DelayedArray_0.34.1         Rhdf5lib_1.30.0            
#>  [33] BiocParallel_1.42.2         jpeg_0.1-11                
#>  [35] data.tree_1.2.0             parallel_4.5.1             
#>  [37] cluster_2.1.8.1             R6_2.6.1                   
#>  [39] bslib_0.9.0                 stringi_1.8.7              
#>  [41] RColorBrewer_1.1-3          GenomicRanges_1.60.0       
#>  [43] jquerylib_0.1.4             SummarizedExperiment_1.38.1
#>  [45] iterators_1.0.14            knitr_1.50                 
#>  [47] IRanges_2.42.0              Matrix_1.7-4               
#>  [49] splines_4.5.1               igraph_2.1.4               
#>  [51] tidyselect_1.2.1            abind_1.4-8                
#>  [53] yaml_2.3.10                 vegan_2.7-1                
#>  [55] codetools_0.2-20            hwriter_1.3.2.1            
#>  [57] lattice_0.22-7              tibble_3.3.0               
#>  [59] plyr_1.8.9                  Biobase_2.68.0             
#>  [61] withr_3.0.2                 ShortRead_1.66.0           
#>  [63] S7_0.2.0                    evaluate_1.0.5             
#>  [65] desc_1.4.3                  survival_3.8-3             
#>  [67] RcppParallel_5.1.11-1       Biostrings_2.76.0          
#>  [69] pillar_1.11.1               MatrixGenerics_1.20.0      
#>  [71] DT_0.34.0                   foreach_1.5.2              
#>  [73] stats4_4.5.1                generics_0.1.4             
#>  [75] S4Vectors_0.46.0            scales_1.4.0               
#>  [77] glue_1.8.0                  tools_4.5.1                
#>  [79] interp_1.1-6                data.table_1.17.8          
#>  [81] GenomicAlignments_1.44.0    fs_1.6.6                   
#>  [83] rhdf5_2.52.1                grid_4.5.1                 
#>  [85] ape_5.8-1                   crosstalk_1.2.2            
#>  [87] latticeExtra_0.6-31         networkD3_0.4.1            
#>  [89] nlme_3.1-168                GenomeInfoDbData_1.2.14    
#>  [91] cli_3.6.5                   textshaping_1.0.3          
#>  [93] S4Arrays_1.8.1              gtable_0.3.6               
#>  [95] sass_0.4.10                 digest_0.6.37              
#>  [97] BiocGenerics_0.54.0         SparseArray_1.8.1          
#>  [99] htmlwidgets_1.6.4           farver_2.1.2               
#> [101] htmltools_0.5.8.1           pkgdown_2.1.3              
#> [103] multtest_2.64.0             lifecycle_1.0.4            
#> [105] httr_1.4.7                  MASS_7.3-65
```
