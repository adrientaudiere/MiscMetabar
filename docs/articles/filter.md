# Filter taxa and samples

``` r
library(MiscMetabar)
library(formattable)
```

### Filter samples

[Phyloseq](https://joey711.github.io/phyloseq/) package already propose
a function to select samples
([`subset_samples()`](https://rdrr.io/pkg/phyloseq/man/subset_samples-methods.html)),
but in some case, the subset internal function is painful. `MiscMetabar`
propose a complementary function (\[subset_samples_pq()\]) which is more
versatile but need to be used with caution because the order of the
condition must match the orders of the samples.

``` r
data(data_fungi)
cond_samp <- grepl("A1", data_fungi@sam_data[["Sample_names"]])
subset_samples_pq(data_fungi, cond_samp)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 9 samples ]
#> sample_data() Sample Data:       [ 9 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```

### Filter Taxa

#### Filter taxa using condition(s)

[Phyloseq](https://joey711.github.io/phyloseq/) package already propose
a function to select samples
([`subset_taxa()`](https://rdrr.io/pkg/phyloseq/man/subset_taxa-methods.html)),
but in some case, the subset internal function is painful. `MiscMetabar`
propose a complementary function (\[subset_taxa_pq()\]) which is more
versatile and is based on the names of the taxa to match the condition
and the taxa in the phyloseq object. In the example code above, we
filter fungi from the Phylum “Basidiomycota” using
[`phyloseq::subset_taxa()`](https://rdrr.io/pkg/phyloseq/man/subset_taxa-methods.html)
and then select only taxa with more than 1000 nb_sequences using
[`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_pq.md).

``` r
df_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
df_basidio_abundant <- subset_taxa_pq(
  df_basidio,
  colSums(df_basidio@otu_table) > 1000
)
```

#### Filter taxa using blast score against a database

In some cases, we want to select only a given clade of taxon. One
solution is to select only taxa with information given by taxonomic
assignment (e.g. using the function
[`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html)).
However, in some cases, this information lead to false positive and
false negative selection. `MiscMetabar` implement an other method
relying on [blast](https://blast.ncbi.nlm.nih.gov/) software (Altschul
*et al.* 1990 & 1997) and database. The idea is to set a cutoff in four
parameters to select only taxa which are close enough to sequences in
the database :

- **id_filter** (default: 90) cut of in identity percent to keep result.
- **bit_score_filter** (default: 50) cut of in bit score to keep result.
  The higher the bit-score, the better the sequence similarity. The
  bit-score is the requires size of a sequence database in which the
  current match could be found just by chance. The bit-score is a log2
  scaled and normalized raw-score. Each increase by one doubles the
  required database size (2bit-score).
- **min_cover_filter** (default: 50) cut of in query cover (%) to keep
  result.
- **e_value_filter** (default: 1e-30) cut of in e-value (%) to keep
  result. The BLAST E-value is the number of expected hits of similar
  quality (score) that could be found just by chance.

``` r
path_db <- system.file("extdata",
  "100_sp_UNITE_sh_general_release_dynamic.fasta",
  package = "MiscMetabar", mustWork = TRUE
)

suppressWarnings(blast_error_or_not <-
  try(system("blastn 2>&1", intern = TRUE), silent = TRUE))

if (!is(blast_error_or_not, "try-error")) {
  df_blast_80 <- filter_asv_blast(df_basidio, fasta_for_db = path_db)
  df_blast_50 <- filter_asv_blast(df_basidio,
    fasta_for_db = path_db,
    id_filter = 50, e_value_filter = 10,
    bit_score_filter = 20, min_cover_filter = 20
  )

  track_formattable <-
    track_wkflow(
      list(
        "raw data" = df_basidio,
        "id_filter = 80" = df_blast_80,
        "id_filter = 50" = df_blast_50
      )
    )
}
```

``` r
if (!is(blast_error_or_not, "try-error")) {
  formattable(
    track_formattable,
    list(
      area(col = nb_sequences) ~ color_bar("cyan", na.rm = TRUE),
      area(col = nb_clusters) ~ normalize_bar("yellowgreen",
        na.rm = TRUE, min = 0.3
      ),
      area(col = nb_samples) ~ normalize_bar("lightpink",
        na.rm = TRUE, min = 0.3
      )
    )
  )
}
```

|                | nb_sequences | nb_clusters | nb_samples |
|:---------------|-------------:|------------:|-----------:|
| raw data       |       784054 |         345 |        185 |
| id_filter = 80 |       249494 |          45 |        127 |
| id_filter = 50 |       780242 |         326 |        172 |

#### Filter taxa using a known taxa for control

To filter out contamination, one solution is to add a proportion of a
known taxa which is not present in the environment of the study. In that
case we can define some threshold for each sample to discard taxon based
on pseudo-abundance. In the example code above, we select taxon using
the ASV_50 as control through 6 different algorithms.

``` r
res_seq <-
  suppressWarnings(
    subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
      method = "cutoff_seq"
    )
  )
res_mixt <-
  suppressWarnings(
    subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
      method = "cutoff_mixt"
    )
  )
res_diff <- suppressWarnings(
  subset_taxa_tax_control(
    data_fungi,
    as.numeric(data_fungi@otu_table[, 50]),
    method = "cutoff_diff",
    min_diff_for_cutoff = 2
  )
)
res_min <-
  suppressWarnings(
    subset_taxa_tax_control(
      data_fungi,
      as.numeric(data_fungi@otu_table[, 50]),
      method = "min",
      min_diff_for_cutoff = 2
    )
  )
res_max <-
  suppressWarnings(
    subset_taxa_tax_control(
      data_fungi,
      as.numeric(data_fungi@otu_table[, 50]),
      method = "max",
      min_diff_for_cutoff = 2
    )
  )
res_mean <-
  suppressWarnings(
    subset_taxa_tax_control(
      data_fungi,
      as.numeric(data_fungi@otu_table[, 50]),
      method = "mean",
      min_diff_for_cutoff = 2
    )
  )
```

``` r
track_formattable <- track_wkflow(list(
  "raw data" = data_fungi,
  "cutoff_seq" = res_seq,
  "cutoff_mixt" = res_mixt,
  "cutoff_diff" = res_diff,
  "min" = res_min,
  "max" = res_max,
  "mean" = res_mean
))

formattable(
  track_formattable,
  list(
    area(col = nb_sequences) ~ color_bar("cyan"),
    area(col = nb_clusters) ~ normalize_bar("yellowgreen",
      na.rm = TRUE,
      min = 0.3
    ),
    area(col = nb_samples) ~ normalize_bar("lightpink", na.rm = TRUE)
  )
)
```

|             | nb_sequences | nb_clusters | nb_samples |
|:------------|-------------:|------------:|-----------:|
| raw data    |      1839124 |        1420 |        185 |
| cutoff_seq  |      1817543 |        1402 |        185 |
| cutoff_mixt |      1774180 |        1279 |        185 |
| cutoff_diff |      1805274 |        1378 |        185 |
| min         |      1836158 |        1417 |        185 |
| max         |      1752217 |        1214 |        185 |
| mean        |      1790569 |        1325 |        185 |

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
#> [1] formattable_0.2.1  MiscMetabar_0.14.4 purrr_1.1.0        dplyr_1.1.4       
#> [5] dada2_1.36.0       Rcpp_1.1.0         ggplot2_4.0.0      phyloseq_1.52.0   
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
#>  [45] knitr_1.50                  mixtools_2.0.0.1           
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
#>  [67] RcppParallel_5.1.11-1       kernlab_0.9-33             
#>  [69] Biostrings_2.76.0           pillar_1.11.1              
#>  [71] MatrixGenerics_1.20.0       foreach_1.5.2              
#>  [73] stats4_4.5.1                plotly_4.11.0              
#>  [75] generics_0.1.4              S4Vectors_0.46.0           
#>  [77] scales_1.4.0                glue_1.8.0                 
#>  [79] lazyeval_0.2.2              tools_4.5.1                
#>  [81] interp_1.1-6                data.table_1.17.8          
#>  [83] GenomicAlignments_1.44.0    fs_1.6.6                   
#>  [85] rhdf5_2.52.1                grid_4.5.1                 
#>  [87] tidyr_1.3.1                 ape_5.8-1                  
#>  [89] latticeExtra_0.6-31         nlme_3.1-168               
#>  [91] GenomeInfoDbData_1.2.14     cli_3.6.5                  
#>  [93] textshaping_1.0.3           segmented_2.1-4            
#>  [95] viridisLite_0.4.2           S4Arrays_1.8.1             
#>  [97] gtable_0.3.6                sass_0.4.10                
#>  [99] digest_0.6.37               BiocGenerics_0.54.0        
#> [101] SparseArray_1.8.1           htmlwidgets_1.6.4          
#> [103] farver_2.1.2                htmltools_0.5.8.1          
#> [105] pkgdown_2.1.3               multtest_2.64.0            
#> [107] lifecycle_1.0.4             httr_1.4.7                 
#> [109] MASS_7.3-65
```
