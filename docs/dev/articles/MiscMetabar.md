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
#>  [7] crosstalk_1.2.2       Rdpack_2.6.6          vctrs_0.7.3          
#> [10] tools_4.6.1           generics_0.1.4        biomformat_1.40.0    
#> [13] stats4_4.6.1          parallel_4.6.1        tibble_3.3.1         
#> [16] cluster_2.1.8.2       pkgconfig_2.0.3       Matrix_1.7-5         
#> [19] data.table_1.18.4     RColorBrewer_1.1-3    S7_0.2.2             
#> [22] desc_1.4.3            S4Vectors_0.50.1      RcppParallel_5.1.11-2
#> [25] lifecycle_1.0.5       compiler_4.6.1        farver_2.1.2         
#> [28] stringr_1.6.0         textshaping_1.0.5     Biostrings_2.80.1    
#> [31] data.tree_1.2.0       Seqinfo_1.2.0         codetools_0.2-20     
#> [34] permute_0.9-10        htmltools_0.5.9       sass_0.4.10          
#> [37] yaml_2.3.12           pkgdown_2.2.0         pillar_1.11.1        
#> [40] crayon_1.5.3          jquerylib_0.1.4       MASS_7.3-65          
#> [43] DT_0.34.0             cachem_1.1.0          vegan_2.7-5          
#> [46] iterators_1.0.14      foreach_1.5.2         nlme_3.1-169         
#> [49] tidyselect_1.2.1      digest_0.6.39         stringi_1.8.7        
#> [52] reshape2_1.4.5        labeling_0.4.3        splines_4.6.1        
#> [55] ade4_1.7-24           fastmap_1.2.0         grid_4.6.1           
#> [58] cli_3.6.6             magrittr_2.0.5        survival_3.8-6       
#> [61] ape_5.8-1             withr_3.0.3           scales_1.4.0         
#> [64] rmarkdown_2.31        XVector_0.52.0        networkD3_0.4.1      
#> [67] igraph_2.3.3          multtest_2.68.0       otel_0.2.0           
#> [70] ragg_1.5.2            evaluate_1.0.5        knitr_1.51           
#> [73] rbibutils_2.4.1       IRanges_2.46.0        mgcv_1.9-4           
#> [76] rlang_1.2.0           divent_0.5-4          Rcpp_1.1.1-1.1       
#> [79] glue_1.8.1            BiocGenerics_0.58.1   jsonlite_2.0.0       
#> [82] R6_2.6.1              plyr_1.8.9            systemfonts_1.3.2    
#> [85] fs_2.1.0
```
