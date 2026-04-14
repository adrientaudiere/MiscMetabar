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
if (packageVersion("ggplot2") < "4.0.0") {
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
#>   [1] bitops_1.0-9                deldir_2.0-4               
#>   [3] permute_0.9-10              rlang_1.1.7                
#>   [5] magrittr_2.0.4              ade4_1.7-23                
#>   [7] otel_0.2.0                  matrixStats_1.5.0          
#>   [9] compiler_4.5.2              mgcv_1.9-4                 
#>  [11] png_0.1-8                   systemfonts_1.3.1          
#>  [13] vctrs_0.7.1                 reshape2_1.4.5             
#>  [15] stringr_1.6.0               pwalign_1.6.0              
#>  [17] pkgconfig_2.0.3             crayon_1.5.3               
#>  [19] fastmap_1.2.0               XVector_0.50.0             
#>  [21] labeling_0.4.3              Rsamtools_2.26.0           
#>  [23] rmarkdown_2.30              ragg_1.5.0                 
#>  [25] xfun_0.56                   cachem_1.1.0               
#>  [27] cigarillo_1.0.0             jsonlite_2.0.0             
#>  [29] biomformat_1.38.0           rhdf5filters_1.22.0        
#>  [31] DelayedArray_0.36.0         Rhdf5lib_1.32.0            
#>  [33] BiocParallel_1.44.0         jpeg_0.1-11                
#>  [35] data.tree_1.2.0             parallel_4.5.2             
#>  [37] cluster_2.1.8.2             R6_2.6.1                   
#>  [39] bslib_0.10.0                stringi_1.8.7              
#>  [41] RColorBrewer_1.1-3          GenomicRanges_1.62.1       
#>  [43] jquerylib_0.1.4             Seqinfo_1.0.0              
#>  [45] SummarizedExperiment_1.40.0 iterators_1.0.14           
#>  [47] knitr_1.51                  IRanges_2.44.0             
#>  [49] Matrix_1.7-4                splines_4.5.2              
#>  [51] igraph_2.2.1                tidyselect_1.2.1           
#>  [53] abind_1.4-8                 yaml_2.3.12                
#>  [55] vegan_2.7-2                 codetools_0.2-20           
#>  [57] hwriter_1.3.2.1             lattice_0.22-9             
#>  [59] tibble_3.3.1                plyr_1.8.9                 
#>  [61] Biobase_2.70.0              withr_3.0.2                
#>  [63] ShortRead_1.68.0            S7_0.2.1                   
#>  [65] evaluate_1.0.5              desc_1.4.3                 
#>  [67] survival_3.8-6              RcppParallel_5.1.11-1      
#>  [69] Biostrings_2.78.0           pillar_1.11.1              
#>  [71] MatrixGenerics_1.22.0       DT_0.34.0                  
#>  [73] foreach_1.5.2               stats4_4.5.2               
#>  [75] generics_0.1.4              S4Vectors_0.48.0           
#>  [77] scales_1.4.0                glue_1.8.0                 
#>  [79] tools_4.5.2                 interp_1.1-6               
#>  [81] data.table_1.18.2.1         GenomicAlignments_1.46.0   
#>  [83] fs_1.6.6                    rhdf5_2.54.1               
#>  [85] grid_4.5.2                  ape_5.8-1                  
#>  [87] crosstalk_1.2.2             latticeExtra_0.6-31        
#>  [89] networkD3_0.4.1             nlme_3.1-168               
#>  [91] cli_3.6.5                   textshaping_1.0.4          
#>  [93] S4Arrays_1.10.1             gtable_0.3.6               
#>  [95] sass_0.4.10                 digest_0.6.39              
#>  [97] BiocGenerics_0.56.0         SparseArray_1.10.8         
#>  [99] htmlwidgets_1.6.4           farver_2.1.2               
#> [101] htmltools_0.5.9             pkgdown_2.2.0              
#> [103] multtest_2.66.0             lifecycle_1.0.5            
#> [105] MASS_7.3-65
```
