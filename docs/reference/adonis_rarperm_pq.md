# Permanova (adonis) on permutations of rarefaction even depth

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Permanova are computed on a given number of rarefaction with different
seed.number. This reduce the risk of a random drawing of a exceptional
situation of an unique rarefaction.

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
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).
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
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html).

- ...:

  Other params to be passed on to
  [`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md)
  function

## Value

A list of three dataframe representing the mean, the minimum quantile
and the maximum quantile value for adonis results. See
[`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md).

## See also

[`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("vegan")) {
  data_fungi_woNA <-
    subset_samples(data_fungi, !is.na(Time) & !is.na(Height))
  adonis_rarperm_pq(data_fungi_woNA, "Time*Height", na_remove = TRUE, nperm = 3)
}
#> 
  |                                                        
  |                                                  |   0%
#> `set.seed(1)` was used to initialize repeatable random subsampling.
#> Please record this for your records so others can reproduce.
#> Try `set.seed(1); .Random.seed` for the full vector
#> ...
#> 1156OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#> Time
#> Height
#> 
  |                                                        
  |=================                                 |  33%
#> `set.seed(2)` was used to initialize repeatable random subsampling.
#> Please record this for your records so others can reproduce.
#> Try `set.seed(2); .Random.seed` for the full vector
#> ...
#> 1160OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#> Time
#> Height
#> 
  |                                                        
  |=================================                 |  67%
#> `set.seed(3)` was used to initialize repeatable random subsampling.
#> Please record this for your records so others can reproduce.
#> Try `set.seed(3); .Random.seed` for the full vector
#> ...
#> 1180OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#> Time
#> Height
#> 
  |                                                        
  |==================================================| 100%
#> $mean
#>           Df  SumOfSqs         R2        F     Pr(>F)
#> Model      5  2.724437 0.05067775 1.121097 0.05733333
#> Residual 105 51.040745 0.94932225       NA         NA
#> Total    110 53.765182 1.00000000       NA         NA
#> 
#> $quantile_min
#>           Df  SumOfSqs        R2        F  Pr(>F)
#> Model      5  2.633518 0.0488029 1.077451 0.03225
#> Residual 105 50.681723 0.9479599       NA      NA
#> Total    110 53.463976 1.0000000       NA      NA
#> 
#> $quantile_max
#>           Df  SumOfSqs         R2        F Pr(>F)
#> Model      5  2.782253 0.05204009 1.152836 0.0997
#> Residual 105 51.328995 0.95119710       NA     NA
#> Total    110 53.962513 1.00000000       NA     NA
#> 
```
