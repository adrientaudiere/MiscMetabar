# Compare samples in pairs using diversity and number of ASV including shared ASV.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

For the moment refseq slot need to be not Null.

## Usage

``` r
compare_pairs_pq(
  physeq = NULL,
  bifactor = NULL,
  modality = NULL,
  merge_sample_by = NULL,
  nb_min_seq = 0,
  veg_index = "shannon",
  na_remove = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- bifactor:

  (required) a factor (present in the `sam_data` slot of the physeq
  object) presenting the pair names

- modality:

  the name of the column in the `sam_data` slot of the physeq object to
  split samples by pairs

- merge_sample_by:

  a vector to determine which samples to merge using the
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  function. Need to be in `physeq@sam_data`

- nb_min_seq:

  minimum number of sequences per sample to count the ASV/OTU

- veg_index:

  (default: `"shannon"`) diversity index. `"shannon"` and `"simpson"`
  are computed via `divent`; other names are forwarded to
  [`vegan::diversity()`](https://vegandevs.github.io/vegan/reference/diversity.html).

- na_remove:

  (logical, default TRUE) If set to TRUE, remove samples with NA in the
  variables set in bifactor, modality and merge_sample_by. NA in
  variables are well managed even if na_remove = FALSE, so na_remove may
  be useless.

- ...:

  Additional arguments passed to
  [`divent::ent_shannon()`](https://ericmarcon.github.io/divent/reference/ent_shannon.html)
  or
  [`divent::ent_simpson()`](https://ericmarcon.github.io/divent/reference/ent_simpson.html)
  when `veg_index` is `"shannon"` or `"simpson"`.

## Value

A tibble with information about the number of shared ASV, shared number
of sequences and diversity

## Examples

``` r
data_fungi_mini_lh <- subset_samples(data_fungi_mini, Height %in% c("Low", "High"))
compare_pairs_pq(data_fungi_mini_lh, bifactor = "Height", merge_sample_by = "Height")
#> Taxa are now in columns.
#> Cleaning suppress 5 taxa and 0 samples.
#> # A tibble: 1 × 13
#>   modality nb_ASV_High nb_ASV_Low nb_shared_ASV div_High div_Low nb_shared_seq
#>   <chr>          <dbl>      <dbl>         <dbl>    <dbl>   <dbl>         <dbl>
#> 1 Height            34         34            28     2.56    2.54        151756
#> # ℹ 6 more variables: percent_shared_seq_High <dbl>,
#> #   percent_shared_seq_Low <dbl>, percent_shared_ASV_High <dbl>,
#> #   percent_shared_ASV_Low <dbl>, ratio_nb_High_Low <dbl>,
#> #   ratio_div_High_Low <dbl>
compare_pairs_pq(data_fungi_mini_lh,
  bifactor = "Height",
  merge_sample_by = "Height", modality = "Time"
)
#> Taxa are now in columns.
#> 7 samples were discarded due to NA in variable modality.
#> Cleaning suppress 5 taxa and 0 samples.
#> # A tibble: 4 × 13
#>   modality nb_ASV_High nb_ASV_Low nb_shared_ASV div_High div_Low nb_shared_seq
#>   <chr>          <dbl>      <dbl>         <dbl>    <dbl>   <dbl>         <dbl>
#> 1 0                 13         17             9     1.6     0.74         22914
#> 2 5                 16         20             7     2.1     1.93         15648
#> 3 10                11         13             8     1.34    1.27         11519
#> 4 15                15         16            10     1.68    1.49         60059
#> # ℹ 6 more variables: percent_shared_seq_High <dbl>,
#> #   percent_shared_seq_Low <dbl>, percent_shared_ASV_High <dbl>,
#> #   percent_shared_ASV_Low <dbl>, ratio_nb_High_Low <dbl>,
#> #   ratio_div_High_Low <dbl>
```
