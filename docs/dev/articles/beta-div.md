# beta-div

``` r
library(MiscMetabar)
library(lefser)
library(ALDEx2)
data(data_fungi)
```

### Permanova

``` r
data_fungi_woNA4height <- subset_samples(data_fungi, !is.na(data_fungi@sam_data$Height))
res_ado <- adonis_pq(data_fungi_woNA4height, "Tree_name+Height")
knitr::kable(res_ado)
```

|          |  Df | SumOfSqs |        R2 |        F | Pr(\>F) |
|:---------|----:|---------:|----------:|---------:|--------:|
| Model    |  63 | 36.92559 | 0.5881754 | 1.518899 |   0.001 |
| Residual |  67 | 25.85431 | 0.4118246 |       NA |      NA |
| Total    | 130 | 62.77990 | 1.0000000 |       NA |      NA |

### Graph Test

``` r
data_fungi_woNA4height <- subset_samples(data_fungi, !is.na(data_fungi@sam_data$Height))
graph_test_pq(data_fungi_woNA4height, "Height")
```

![](beta-div_files/figure-html/unnamed-chunk-3-1.png)

### Circle of ASVs

``` r
circle_pq(data_fungi_woNA4height, "Height")
```

![](beta-div_files/figure-html/unnamed-chunk-4-1.png)

### Ordination

#### PCoA

``` r
plot_ordination_pq(data_fungi, method = "robust.aitchison", ordination_method = "PCoA", color = "Height")
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the phyloseq package.
#>   Please report the issue at <https://github.com/joey711/phyloseq/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

![](beta-div_files/figure-html/unnamed-chunk-5-1.png)

``` r
plot_ordination(data_fungi,
  ordination =
    ordinate(data_fungi, method = "PCoA", distance = "bray"), color = "Height"
)
```

![](beta-div_files/figure-html/unnamed-chunk-5-2.png)

#### NMDS

``` r
plot_ordination_pq(data_fungi, method = "robust.aitchison", color = "Height") +
  plot_ordination_pq(data_fungi, method = "bray", color = "Height")
#> Run 0 stress 0.1133285 
#> Run 1 stress 0.2031998 
#> Run 2 stress 0.1425415 
#> Run 3 stress 0.1273063 
#> Run 4 stress 0.1950973 
#> Run 5 stress 0.1620808 
#> Run 6 stress 0.1595056 
#> Run 7 stress 0.1796361 
#> Run 8 stress 0.1721042 
#> Run 9 stress 0.1980584 
#> Run 10 stress 0.1284201 
#> Run 11 stress 0.1971293 
#> Run 12 stress 0.2005546 
#> Run 13 stress 0.168934 
#> Run 14 stress 0.1812448 
#> Run 15 stress 0.1881208 
#> Run 16 stress 0.1575817 
#> Run 17 stress 0.162739 
#> Run 18 stress 0.132636 
#> Run 19 stress 0.152047 
#> Run 20 stress 0.2041781 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     13: stress ratio > sratmax
#>      7: scale factor of the gradient < sfgrmin
#> Run 0 stress 0.2437839 
#> Run 1 stress 0.2485997 
#> Run 2 stress 0.2472849 
#> Run 3 stress 0.2450257 
#> Run 4 stress 0.2394284 
#> ... New best solution
#> ... Procrustes: rmse 0.04540589  max resid 0.1926736 
#> Run 5 stress 0.2471695 
#> Run 6 stress 0.2513681 
#> Run 7 stress 0.2429278 
#> Run 8 stress 0.2489632 
#> Run 9 stress 0.2489441 
#> Run 10 stress 0.2453355 
#> Run 11 stress 0.2466278 
#> Run 12 stress 0.2455799 
#> Run 13 stress 0.2483674 
#> Run 14 stress 0.2460726 
#> Run 15 stress 0.2460076 
#> Run 16 stress 0.246254 
#> Run 17 stress 0.2453415 
#> Run 18 stress 0.2471357 
#> Run 19 stress 0.2432681 
#> Run 20 stress 0.2522929 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>      3: no. of iterations >= maxit
#>     17: stress ratio > sratmax
```

![](beta-div_files/figure-html/unnamed-chunk-6-1.png)

#### TSNE

``` r
plot_tsne_pq(data_fungi, fact = "Height")
```

![](beta-div_files/figure-html/unnamed-chunk-7-1.png)

#### UMAP

``` r
df_umap <- umap_pq(data_fungi)
ggplot(df_umap, aes(x = x_umap, y = y_umap, col = Height)) +
  geom_point(size = 2)
```

![](beta-div_files/figure-html/unnamed-chunk-8-1.png)

### Compare two (group of) samples

#### Biplot

``` r
data_fungi_low_high <- subset_samples(
  data_fungi,
  data_fungi@sam_data$Height %in%
    c("Low", "High")
)
data_fungi_low_high <- subset_taxa_pq(
  data_fungi_low_high,
  taxa_sums(data_fungi_low_high) > 5000
)
biplot_pq(data_fungi_low_high, fact = "Height", merge_sample_by = "Height")
```

![](beta-div_files/figure-html/unnamed-chunk-9-1.png)

#### Compare two (group of) samples with a table

``` r
compare_pairs_pq(data_fungi_low_high,
  bifactor = "Height",
  merge_sample_by = "Height",
  modality = "Time"
)
#> # A tibble: 4 × 13
#>   modality nb_ASV_High nb_ASV_Low nb_shared_ASV div_High div_Low nb_shared_seq
#>   <chr>          <dbl>      <dbl>         <dbl>    <dbl>   <dbl>         <dbl>
#> 1 0                 12         16             9     1.8     1.37         57639
#> 2 5                 20         18            14     1.95    1.98         76006
#> 3 10                11         13            10     1.18    1.25         47042
#> 4 15                17         19            12     2       2.04        161348
#> # ℹ 6 more variables: percent_shared_seq_High <dbl>,
#> #   percent_shared_seq_Low <dbl>, percent_shared_ASV_High <dbl>,
#> #   percent_shared_ASV_Low <dbl>, ratio_nb_High_Low <dbl>,
#> #   ratio_div_High_Low <dbl>
```

### Venn diagram

``` r
library("grid")
venn_pq(data_fungi, fact = "Height")
```

![](beta-div_files/figure-html/unnamed-chunk-11-1.png)

``` r
ggvenn_pq(data_fungi, fact = "Height") +
  ggplot2::scale_fill_distiller(palette = "BuPu", direction = 1) +
  labs(title = "Share number of ASV among Height in tree")
```

![](beta-div_files/figure-html/unnamed-chunk-12-1.png)

``` r
ggvenn_pq(data_fungi, fact = "Height", min_nb_seq = 5000) +
  ggplot2::scale_fill_distiller(palette = "BuPu", direction = 1) +
  labs(title = "Share number of ASV with more than 5000 seqs")
```

![](beta-div_files/figure-html/unnamed-chunk-13-1.png)

``` r
ggvenn_pq(data_fungi,
  fact = "Height", taxonomic_rank = "Genus",
  min_nb_seq = 100
) +
  ggplot2::scale_fill_distiller(palette = "BuPu", direction = 1) +
  labs(title = "Share number of Genus represented by at least one ASV with more than 100 seqs")
```

![](beta-div_files/figure-html/unnamed-chunk-14-1.png)

### Upset plot

Venn diagram can quickly become complex to read when the number of
modalities increase. One graphical solution is upset plot. MiscMetabar
propose a solution based on the package
[ComplexUpset](https://krassowski.github.io/complex-upset/).

``` r
upset_pq(data_fungi, fact = "Height")
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ComplexUpset package.
#>   Please report the issue at
#>   <https://github.com/krassowski/complex-upset/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

![](beta-div_files/figure-html/unnamed-chunk-15-1.png)

``` r
upset_pq(data_fungi, fact = "Time")
```

![](beta-div_files/figure-html/unnamed-chunk-16-1.png)

`ComplexUpset` package allow powerful configuration of you plot as you
can see in the following figure.

``` r
upset_pq(
  data_fungi,
  fact = "Time",
  width_ratio = 0.2,
  annotations = list(
    "Sequences per ASV \n (log10)" = (
      ggplot(mapping = aes(y = log10(Abundance)))
      +
        geom_jitter(aes(
          color =
            Abundance
        ), na.rm = TRUE)
        +
        geom_violin(alpha = 0.5, na.rm = TRUE) +
        theme(legend.key.size = unit(0.2, "cm")) +
        theme(axis.text = element_text(size = 12))
    ),
    "ASV per phylum" = (
      ggplot(mapping = aes(fill = Phylum))
      +
        geom_bar() +
        ylab("ASV per phylum") +
        theme(legend.key.size = unit(0.2, "cm")) +
        theme(axis.text = element_text(size = 12))
    )
  )
)
```

![](beta-div_files/figure-html/unnamed-chunk-17-1.png)

### Change in abundance across a factor

#### Benchdamic

There is a lot of available methods. Please refer to R package
[benchdamic](https://github.com/mcalgaro93/benchdamic) for a list of
method and a implementation of a benchmark for your data.

See also [Gamboa-Tuz et
al.](https://www.biorxiv.org/content/biorxiv/early/2025/02/17/2025.02.13.638109.full.pdf)
who claimed that “\[…\] compositional DA methods are not beneficial but
rather lack sensitivity, show increased variability in
constant-abundance spike-ins, and, most surprisingly, more frequently
produce paradoxical results with DA in the wrong direction for the
low-diversity microbiome. Conversely, commonly used methods in
microbiome literature, such as LEfSe, the Wilcoxon test, and
RNA-seq-derived methods, performed best.”

#### Library requirement for Debian Linux OS

``` sh
sudo apt-get install libgsl-dev libmpfr-dev
```

#### Using Deseq2 package

``` r
data("GlobalPatterns", package = "phyloseq")
GP <- subset_samples(
  GlobalPatterns,
  GlobalPatterns@sam_data$SampleType %in% c("Soil", "Skin")
)

plot_deseq2_pq(GP, c("SampleType", "Soil", "Skin"), pval = 0.001)
```

![](beta-div_files/figure-html/unnamed-chunk-18-1.png)

#### Using Linear discriminant analysis (LDA) Effect Size (LEfSe)

``` r
res_lefse <- lefser_pq(data_fungi, bifactor = "Height", modalities = c("Low", "High"))
lefser::lefserPlot(res_lefse)
```

![](beta-div_files/figure-html/unnamed-chunk-19-1.png)

``` r

res_lefse_clade <- lefser_pq(data_fungi, bifactor = "Height", modalities = c("Low", "High"), by_clade = TRUE)
lefser::lefserPlot(res_lefse_clade)
```

![](beta-div_files/figure-html/unnamed-chunk-19-2.png)

#### Using ALDEx2 package

``` r
res_aldex <- aldex_pq(data_fungi_mini, bifactor = "Height", modalities = c("Low", "High"))

ALDEx2::aldex.plot(res_aldex, type = "volcano")
```

![](beta-div_files/figure-html/unnamed-chunk-20-1.png)

#### Using ancombc

``` r
res_ancombc <- ancombc_pq(
  data_fungi_mini,
  fact = "Height",
  levels_fact = c("Low", "High"),
  verbose = TRUE
)
#> Warning: The group variable has < 3 categories 
#> The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
#> Warning: The number of taxa used for estimating sample-specific biases is: 11
#> A large number of taxa (>50) is required for the consistent estimation of biases
#> Warning: Estimation of sampling fractions failed for the following samples:
#> B18-006-B_S19_MERGED.fastq.gz, DY5-004-B_S96_MERGED.fastq.gz, W26-001-B_S165_MERGED.fastq.gz, X24-009-B_S170_MERGED.fastq.gz, X29-004-B_S174_MERGED.fastq.gz, Y28-002-B_S178_MERGED.fastq.gz, Z29-001-H_S185_MERGED.fastq.gz
#> These samples may have an excessive number of zero values
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
#>  [1] LC_CTYPE=en_US.UTF-8          LC_NUMERIC=C                 
#>  [3] LC_TIME=en_US.UTF-8           LC_COLLATE=en_US.UTF-8       
#>  [5] LC_MONETARY=en_US.UTF-8       LC_MESSAGES=en_US.UTF-8      
#>  [7] LC_PAPER=en_US.UTF-8          LC_NAME=en_US.UTF-8          
#>  [9] LC_ADDRESS=en_US.UTF-8        LC_TELEPHONE=en_US.UTF-8     
#> [11] LC_MEASUREMENT=en_US.UTF-8    LC_IDENTIFICATION=en_US.UTF-8
#> 
#> time zone: Europe/Paris
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] grid      stats4    stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] doRNG_1.8.6.3               rngtools_1.5.2             
#>  [3] foreach_1.5.2               ALDEx2_1.44.0              
#>  [5] latticeExtra_0.6-31         lattice_0.22-9             
#>  [7] zCompositions_1.6.1         survival_3.8-6             
#>  [9] truncnorm_1.0-9             MASS_7.3-65                
#> [11] lefser_1.22.0               SummarizedExperiment_1.42.0
#> [13] Biobase_2.72.0              GenomicRanges_1.64.0       
#> [15] Seqinfo_1.2.0               IRanges_2.46.0             
#> [17] S4Vectors_0.50.1            BiocGenerics_0.58.1        
#> [19] generics_0.1.4              MatrixGenerics_1.24.0      
#> [21] matrixStats_1.5.0           MiscMetabar_0.17.0.9000    
#> [23] dplyr_1.2.1                 ggplot2_4.0.3              
#> [25] phyloseq_1.56.0            
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_2.1.0                        DirichletMultinomial_1.54.0    
#>   [3] doParallel_1.0.17               httr_1.4.8                     
#>   [5] RColorBrewer_1.1-3              numDeriv_2016.8-1.1            
#>   [7] backports_1.5.1                 tools_4.6.1                    
#>   [9] utf8_1.2.6                      R6_2.6.1                       
#>  [11] vegan_2.7-5                     lazyeval_0.2.3                 
#>  [13] mgcv_1.9-4                      permute_0.9-10                 
#>  [15] withr_3.0.3                     gridExtra_2.3                  
#>  [17] cli_3.6.6                       textshaping_1.0.5              
#>  [19] network_1.20.0                  sandwich_3.1-1                 
#>  [21] labeling_0.4.3                  sass_0.4.10                    
#>  [23] mvtnorm_1.4-1                   S7_0.2.2                       
#>  [25] readr_2.2.0                     proxy_0.4-29                   
#>  [27] askpass_1.2.1                   pkgdown_2.2.0                  
#>  [29] systemfonts_1.3.2               yulab.utils_0.2.4              
#>  [31] foreign_0.8-91                  scater_1.40.1                  
#>  [33] decontam_1.32.0                 readxl_1.5.0                   
#>  [35] rstudioapi_0.19.0               gridGraphics_0.5-1             
#>  [37] ggVennDiagram_1.5.7             shape_1.4.6.1                  
#>  [39] gtools_3.9.5                    Matrix_1.7-5                   
#>  [41] interp_1.1-6                    biomformat_1.40.0              
#>  [43] ggbeeswarm_0.7.3                DescTools_0.99.60              
#>  [45] DECIPHER_3.8.0                  abind_1.4-8                    
#>  [47] lifecycle_1.0.5                 multcomp_1.4-30                
#>  [49] yaml_2.3.12                     SparseArray_1.12.2             
#>  [51] Rtsne_0.17                      crayon_1.5.3                   
#>  [53] haven_2.5.5                     beachmat_2.28.0                
#>  [55] ComplexUpset_1.3.3              sna_2.8                        
#>  [57] venneuler_1.1-4                 pillar_1.11.1                  
#>  [59] knitr_1.51                      boot_1.3-32                    
#>  [61] gld_2.6.8                       codetools_0.2-20               
#>  [63] glue_1.8.1                      ggiraph_0.9.6                  
#>  [65] ggfun_0.2.0                     fontLiberation_0.1.0           
#>  [67] data.table_1.18.4               MultiAssayExperiment_1.38.0    
#>  [69] vctrs_0.7.3                     png_0.1-9                      
#>  [71] treeio_1.36.1                   Rdpack_2.6.6                   
#>  [73] cellranger_1.1.0                testthat_3.3.2                 
#>  [75] gtable_0.3.6                    cachem_1.1.0                   
#>  [77] zigg_0.0.2                      xfun_0.58                      
#>  [79] rbibutils_2.4.1                 S4Arrays_1.12.0                
#>  [81] Rfast_2.1.5.2                   libcoin_1.0-12                 
#>  [83] reformulas_0.4.4                coda_0.19-4.1                  
#>  [85] SingleCellExperiment_1.34.0     rJava_1.0-18                   
#>  [87] iterators_1.0.14                bluster_1.22.0                 
#>  [89] directlabels_2026.4.23          TH.data_1.1-5                  
#>  [91] nlme_3.1-169                    ANCOMBC_2.14.0                 
#>  [93] phyloseqGraphTest_0.1.1         ggtree_4.2.0                   
#>  [95] fontquiver_0.2.1                bslib_0.11.0                   
#>  [97] irlba_2.3.7                     rpart_4.1.27                   
#>  [99] vipor_0.4.7                     otel_0.2.0                     
#> [101] Hmisc_5.2-5                     colorspace_2.1-2               
#> [103] DBI_1.3.0                       nnet_7.3-20                    
#> [105] ade4_1.7-24                     Exact_3.3                      
#> [107] DESeq2_1.52.0                   tidyselect_1.2.1               
#> [109] compiler_4.6.1                  microbiome_1.34.0              
#> [111] htmlTable_2.5.0                 BiocNeighbors_2.6.0            
#> [113] expm_1.0-0                      desc_1.4.3                     
#> [115] fontBitstreamVera_0.1.1         DelayedArray_0.38.2            
#> [117] checkmate_2.3.4                 scales_1.4.0                   
#> [119] quadprog_1.5-8                  rappdirs_0.3.4                 
#> [121] stringr_1.6.0                   digest_0.6.39                  
#> [123] minqa_1.2.8                     rmarkdown_2.31                 
#> [125] XVector_0.52.0                  base64enc_0.1-6                
#> [127] htmltools_0.5.9                 pkgconfig_2.0.3                
#> [129] jpeg_0.1-11                     lme4_2.0-1                     
#> [131] umap_0.2.10.0                   sparseMatrixStats_1.24.0       
#> [133] fastmap_1.2.0                   rlang_1.2.0                    
#> [135] GlobalOptions_0.1.4             htmlwidgets_1.6.4              
#> [137] DelayedMatrixStats_1.34.0       farver_2.1.2                   
#> [139] jquerylib_0.1.4                 energy_1.7-12                  
#> [141] zoo_1.8-15                      jsonlite_2.0.0                 
#> [143] BiocParallel_1.46.0             statnet.common_4.13.0          
#> [145] BiocSingular_1.28.0             magrittr_2.0.5                 
#> [147] Formula_1.2-5                   modeltools_0.2-24              
#> [149] scuttle_1.22.0                  ggnetwork_0.5.14               
#> [151] ggplotify_0.1.3                 patchwork_1.3.2                
#> [153] Rcpp_1.1.1-1.1                  ape_5.8-1                      
#> [155] viridis_0.6.5                   gdtools_0.5.1                  
#> [157] reticulate_1.46.0               stringi_1.8.7                  
#> [159] rootSolve_1.8.2.4               brio_1.1.5                     
#> [161] plyr_1.8.9                      parallel_4.6.1                 
#> [163] ggrepel_0.9.8                   forcats_1.0.1                  
#> [165] lmom_3.3                        deldir_2.0-4                   
#> [167] Biostrings_2.80.1               splines_4.6.1                  
#> [169] multtest_2.68.0                 hms_1.1.4                      
#> [171] circlize_0.4.18                 locfit_1.5-9.12                
#> [173] igraph_2.3.3                    reshape2_1.4.5                 
#> [175] ScaledMatrix_1.20.0             evaluate_1.0.5                 
#> [177] RcppParallel_5.1.11-2           nloptr_2.2.1                   
#> [179] tzdb_0.5.0                      tidyr_1.3.2                    
#> [181] openssl_2.4.2                   purrr_1.2.2                    
#> [183] ecodive_2.2.6                   rsvd_1.0.5                     
#> [185] coin_1.4-3                      divent_0.5-4                   
#> [187] e1071_1.7-17                    RSpectra_0.16-2                
#> [189] tidytree_0.4.7                  class_7.3-23                   
#> [191] viridisLite_0.4.3               ragg_1.5.2                     
#> [193] gsl_2.1-9                       tibble_3.3.1                   
#> [195] lmerTest_3.2-1                  aplot_0.2.9                    
#> [197] beeswarm_0.4.0                  cluster_2.1.8.2                
#> [199] TreeSummarizedExperiment_2.20.0 mia_1.20.0
```
