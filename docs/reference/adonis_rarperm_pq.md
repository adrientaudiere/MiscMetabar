# Permanova (adonis) on permutations of rarefaction even depth

[![lifecycle-deprecated](https://img.shields.io/badge/lifecycle-deprecated-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Usage

``` r
adonis_rarperm_pq(
  physeq,
  formula,
  dist_method = "bray",
  merge_sample_by = NULL,
  na_remove = FALSE,
  rarefy_nb_seqs = FALSE,
  verbose = TRUE,
  nperm = 99,
  progress_bar = TRUE,
  quantile_prob = 0.975,
  sample.size = min(sample_sums(physeq)),
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- formula:

  (required) the right part of a formula for
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html).
  Variables must be present in the `physeq@sam_data` slot.

- dist_method:

  (default "bray") the distance used. See
  [`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
  for all available distances or run
  [`phyloseq::distanceMethodList()`](https://rdrr.io/pkg/phyloseq/man/distanceMethodList.html).
  For aitchison and robust.aitchison distance,
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function is directly used.

- merge_sample_by:

  a vector to determine which samples to merge using the
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  function. Need to be in `physeq@sam_data`

- na_remove:

  (logical, default FALSE) If set to TRUE, remove samples with NA in the
  variables set in formula.

- rarefy_nb_seqs:

  (logical, default FALSE) Rarefy each sample (before merging if
  merge_sample_by is set) using
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html).
  if `correction_for_sample_size` is TRUE, rarefy_nb_seqs will have no
  effect.

- verbose:

  (logical, default TRUE) If TRUE, prompt some messages.

- nperm:

  (int, default = 99) The number of permutations to perform.

- progress_bar:

  (logical, default TRUE) Do we print progress during the calculation.

- quantile_prob:

  (float, `[0:1]`) the value to compute the quantile. Minimum quantile
  is computed using 1-quantile_prob.

- sample.size:

  (int) A single integer value equal to the number of reads being
  simulated, also known as the depth. See
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  and
  [`rarefy_even_depth_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rarefy_even_depth_pq.md).

- ...:

  Other params to be passed on to
  [`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md)
  function

## Value

A list of three dataframe representing the mean, the minimum quantile
and the maximum quantile value for adonis results. See
[`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md).

## NGNjpel6HZEIUCifveBCKjS4TkiX0cqp-1- This function has moved to the `bootpq`

## package. Please use `bootpq::adonis_rarperm_pq()` instead.

Permanova are computed on a given number of rarefaction with different
seed.number. This reduce the risk of a random drawing of a exceptional
situation of an unique rarefaction.

## See also

[`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("vegan")) {
  data_fungi_woNA <-
    subset_samples(data_fungi_mini, !is.na(Time) & !is.na(Height))
  adonis_rarperm_pq(data_fungi_woNA, "Time*Height", na_remove = TRUE, nperm = 3)
}
#>   |                                                          |                                                  |   0%
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#>   |                                                          |=================                                 |  33%
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#>   |                                                          |=================================                 |  67%
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#>   |                                                          |==================================================| 100%
#> $mean
#>          Df  SumOfSqs         R2         F    Pr(>F)
#> Model     5  2.346505 0.06649319 0.9844575 0.5556667
#> Residual 69 32.942384 0.93350681        NA        NA
#> Total    74 35.288889 1.00000000        NA        NA
#> 
#> $quantile_min
#>          Df  SumOfSqs         R2         F Pr(>F)
#> Model     5  2.070369 0.05846765 0.8569665 0.1065
#> Residual 69 32.531964 0.92132923        NA     NA
#> Total    74 35.130000 1.00000000        NA     NA
#> 
#> $quantile_max
#>          Df  SumOfSqs         R2        F  Pr(>F)
#> Model     5  2.778036 0.07867077 1.178625 0.86175
#> Residual 69 33.340964 0.94153235       NA      NA
#> Total    74 35.421333 1.00000000       NA      NA
#> 
# }
```
