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
[`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/subset_taxa_pq.md).

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
#> [1] formattable_0.2.1.9000  MiscMetabar_0.16.1.9000 divent_0.5-3           
#> [4] purrr_1.2.2             dplyr_1.2.1             dada2_1.38.0           
#> [7] Rcpp_1.1.1              ggplot2_4.0.2           phyloseq_1.54.2        
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
#>  [47] knitr_1.51                  mixtools_2.0.0.1           
#>  [49] IRanges_2.44.0              Matrix_1.7-4               
#>  [51] splines_4.5.2               igraph_2.2.3               
#>  [53] tidyselect_1.2.1            abind_1.4-8                
#>  [55] yaml_2.3.12                 vegan_2.7-3                
#>  [57] codetools_0.2-20            hwriter_1.3.2.1            
#>  [59] lattice_0.22-9              tibble_3.3.1               
#>  [61] plyr_1.8.9                  Biobase_2.70.0             
#>  [63] withr_3.0.2                 ShortRead_1.68.0           
#>  [65] S7_0.2.1                    evaluate_1.0.5             
#>  [67] desc_1.4.3                  survival_3.8-6             
#>  [69] RcppParallel_5.1.11-2       kernlab_0.9-33             
#>  [71] Biostrings_2.78.0           pillar_1.11.1              
#>  [73] MatrixGenerics_1.22.0       foreach_1.5.2              
#>  [75] stats4_4.5.2                plotly_4.12.0              
#>  [77] generics_0.1.4              S4Vectors_0.48.1           
#>  [79] scales_1.4.0                glue_1.8.0                 
#>  [81] lazyeval_0.2.3              tools_4.5.2                
#>  [83] interp_1.1-6                data.table_1.18.2.1        
#>  [85] GenomicAlignments_1.46.0    fs_2.0.1                   
#>  [87] rhdf5_2.54.1                grid_4.5.2                 
#>  [89] tidyr_1.3.2                 ape_5.8-1                  
#>  [91] rbibutils_2.4.1             latticeExtra_0.6-31        
#>  [93] nlme_3.1-168                cli_3.6.6                  
#>  [95] textshaping_1.0.5           segmented_2.2-1            
#>  [97] viridisLite_0.4.3           S4Arrays_1.10.1            
#>  [99] gtable_0.3.6                sass_0.4.10                
#> [101] digest_0.6.39               BiocGenerics_0.56.0        
#> [103] SparseArray_1.10.10         htmlwidgets_1.6.4          
#> [105] farver_2.1.2                htmltools_0.5.9            
#> [107] pkgdown_2.2.0               multtest_2.66.0            
#> [109] lifecycle_1.0.5             httr_1.4.8                 
#> [111] MASS_7.3-65
```
