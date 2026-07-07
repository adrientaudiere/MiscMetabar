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
#> [1] formattable_0.2.1       MiscMetabar_0.17.0.9000 dplyr_1.2.1            
#> [4] ggplot2_4.0.3           phyloseq_1.56.0        
#> 
#> loaded via a namespace (and not attached):
#>  [1] ade4_1.7-24           tidyselect_1.2.1      viridisLite_0.4.3    
#>  [4] farver_2.1.2          Biostrings_2.80.1     mixtools_2.0.0.1     
#>  [7] S7_0.2.2              divent_0.5-4          fastmap_1.2.0        
#> [10] lazyeval_0.2.3        digest_0.6.39         lifecycle_1.0.5      
#> [13] cluster_2.1.8.2       survival_3.8-6        kernlab_0.9-33       
#> [16] magrittr_2.0.5        compiler_4.6.1        rlang_1.2.0          
#> [19] sass_0.4.10           tools_4.6.1           igraph_2.3.3         
#> [22] yaml_2.3.12           data.table_1.18.4     knitr_1.51           
#> [25] htmlwidgets_1.6.4     plyr_1.8.9            RColorBrewer_1.1-3   
#> [28] withr_3.0.3           purrr_1.2.2           BiocGenerics_0.58.1  
#> [31] desc_1.4.3            grid_4.6.1            stats4_4.6.1         
#> [34] multtest_2.68.0       biomformat_1.40.0     scales_1.4.0         
#> [37] iterators_1.0.14      MASS_7.3-65           cli_3.6.6            
#> [40] rmarkdown_2.31        vegan_2.7-5           crayon_1.5.3         
#> [43] ragg_1.5.2            generics_0.1.4        otel_0.2.0           
#> [46] RcppParallel_5.1.11-2 httr_1.4.8            reshape2_1.4.5       
#> [49] pbapply_1.7-4         ape_5.8-1             cachem_1.1.0         
#> [52] stringr_1.6.0         splines_4.6.1         parallel_4.6.1       
#> [55] XVector_0.52.0        vctrs_0.7.3           Matrix_1.7-5         
#> [58] jsonlite_2.0.0        IRanges_2.46.0        S4Vectors_0.50.1     
#> [61] systemfonts_1.3.2     foreach_1.5.2         plotly_4.12.0        
#> [64] tidyr_1.3.2           jquerylib_0.1.4       glue_1.8.1           
#> [67] pkgdown_2.2.0         codetools_0.2-20      stringi_1.8.7        
#> [70] gtable_0.3.6          tibble_3.3.1          pillar_1.11.1        
#> [73] htmltools_0.5.9       Seqinfo_1.2.0         R6_2.6.1             
#> [76] textshaping_1.0.5     Rdpack_2.6.6          evaluate_1.0.5       
#> [79] lattice_0.22-9        Biobase_2.72.0        rbibutils_2.4.1      
#> [82] segmented_2.2-1       bslib_0.11.0          Rcpp_1.1.1-1.1       
#> [85] nlme_3.1-169          permute_0.9-10        mgcv_1.9-4           
#> [88] xfun_0.58             fs_2.1.0              pkgconfig_2.0.3
```
