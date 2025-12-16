# Partition the Variation of a phyloseq object with rarefaction permutations

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This is an extension of the function
[`var_par_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_pq.md).
The main addition is the computation of nperm permutations with
rarefaction even depth by sample. The return object

## Usage

``` r
var_par_rarperm_pq(
  physeq,
  list_component,
  dist_method = "bray",
  nperm = 99,
  quantile_prob = 0.975,
  dbrda_computation = FALSE,
  dbrda_signif_pval = 0.05,
  sample.size = min(sample_sums(physeq)),
  verbose = FALSE,
  progress_bar = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- list_component:

  (required) A named list of 2, 3 or four vectors with names from the
  `@sam_data` slot.

- dist_method:

  (default "bray") the distance used. See
  [`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
  for all available distances or run
  [`phyloseq::distanceMethodList()`](https://rdrr.io/pkg/phyloseq/man/distanceMethodList.html).
  For aitchison and robust.aitchison distance,
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function is directly used.#' @param fill_bg

- nperm:

  (int) The number of permutations to perform.

- quantile_prob:

  (float, `[0:1]`) the value to compute the quantile. Minimum quantile
  is compute using 1-quantile_prob.

- dbrda_computation:

  (logical) Do dbrda computations are runned for each individual
  component (each name of the list component) ?

- dbrda_signif_pval:

  (float, `[0:1]`) The value under which the dbrda is considered
  significant.

- sample.size:

  (int) A single integer value equal to the number of reads being
  simulated, also known as the depth. See
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html).

- verbose:

  (logical). If TRUE, print additional information.

- progress_bar:

  (logical, default TRUE) Do we print progress during the calculation?

## Value

A list of class varpart with additional information in the
`$part$indfract` part. Adj.R.square is the mean across permutation.
Adj.R.squared_quantil_min and Adj.R.squared_quantil_max represent the
quantile values of adjusted R squared

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html)
if you use this function.

## See also

[`var_par_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_pq.md),
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html),
[`plot_var_part_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_var_part_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("vegan")) {
  data_fungi_woNA <- subset_samples(data_fungi, !is.na(Time) & !is.na(Height))
  res_var_9 <- var_par_rarperm_pq(
    data_fungi_woNA,
    list_component = list(
      "Time" = c("Time"),
      "Size" = c("Height", "Diameter")
    ),
    nperm = 9,
    dbrda_computation = TRUE
  )
  res_var_2 <- var_par_rarperm_pq(
    data_fungi_woNA,
    list_component = list(
      "Time" = c("Time"),
      "Size" = c("Height", "Diameter")
    ),
    nperm = 2,
    dbrda_computation = TRUE
  )
}
#> 
  |                                                        
  |                                                  |   0%
  |                                                        
  |======                                            |  11%
  |                                                        
  |===========                                       |  22%
  |                                                        
  |=================                                 |  33%
  |                                                        
  |======================                            |  44%
  |                                                        
  |============================                      |  56%
  |                                                        
  |=================================                 |  67%
  |                                                        
  |=======================================           |  78%
  |                                                        
  |============================================      |  89%
  |                                                        
  |==================================================| 100%
  |                                                        
  |                                                  |   0%
  |                                                        
  |=========================                         |  50%
  |                                                        
  |==================================================| 100%
# }
```
