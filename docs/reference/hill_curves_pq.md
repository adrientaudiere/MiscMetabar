# Hill Diversities and Corresponding Accumulation Curves for phyloseq

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wrapper of
[`vegan::renyi()`](https://vegandevs.github.io/vegan/reference/renyi.html)
and
[`vegan::renyiaccum()`](https://vegandevs.github.io/vegan/reference/renyi.html)
functions

## Usage

``` r
hill_curves_pq(
  physeq,
  merge_sample_by = NULL,
  color_fac = NULL,
  hill_scales = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf),
  nperm = NULL,
  na_remove = TRUE,
  wrap_factor = TRUE,
  plot_legend = TRUE,
  linewidth = 2,
  size_point = 2,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- merge_sample_by:

  a vector to determine which samples to merge using the
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  function. Need to be in `physeq@sam_data`

- color_fac:

  (optional): The variable to color the barplot. For ex. same as fact.
  If merge_sample_by is set, color_fac must be nested in the
  merge_sample_by factor. See examples.

- hill_scales:

  Scales of Rényi diversity.

- nperm:

  (int Default NULL) If a integer is set to nperm, nperm permutation are
  computed to draw confidence interval for each curves. The function use
  [`vegan::renyi()`](https://vegandevs.github.io/vegan/reference/renyi.html)
  if nperm is NULL and
  [`vegan::renyiaccum()`](https://vegandevs.github.io/vegan/reference/renyi.html)
  else.

- na_remove:

  (logical, default FALSE) If set to TRUE, remove samples with NA in the
  variables set in merge_sample_by. Not used if merge_sample_by is NULL.

- wrap_factor:

  (logical, default TRUE) Do the plot is wrap by the factor

- plot_legend:

  (logical, default TRUE) If set to FALSE, no legend are plotted.

- linewidth:

  (int, default 2) The linewidth of lines.

- size_point:

  (int, default 1) The size of the point.

- ...:

  Additional arguments passed on to
  [`vegan::renyi()`](https://vegandevs.github.io/vegan/reference/renyi.html)
  function or
  [`vegan::renyiaccum()`](https://vegandevs.github.io/vegan/reference/renyi.html)
  if nperm is not NULL.

## Value

A ggplot2 object

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`vegan::renyi()`](https://vegandevs.github.io/vegan/reference/renyi.html)
or
[`vegan::renyiaccum()`](https://vegandevs.github.io/vegan/reference/renyi.html)
functions

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("vegan")) {
  hill_curves_pq(data_fungi_mini, merge_sample_by = "Time")
  hill_curves_pq(data_fungi_mini, color_fac = "Time", plot_legend = FALSE)
  hill_curves_pq(data_fungi_mini,
    color_fac = "Time", plot_legend = FALSE,
    nperm = 9, size_point = 1, linewidth = 0.5
  )

  hill_curves_pq(data_fungi_mini,
    nperm = 9, plot_legend = FALSE, size_point = 1,
    linewidth = 0.5
  )
  hill_curves_pq(data_fungi_mini, "Height",
    hill_scales = c(0, 1, 2, 8), plot_legend = FALSE
  )
  hill_curves_pq(data_fungi_mini, "Height",
    hill_scales = c(0, 0.5, 1, 2, 4, 8),
    nperm = 9
  )
  hill_curves_pq(data_fungi_mini, "Height", nperm = 9, wrap_factor = FALSE)

  data_fungi_mini@sam_data$H_T <- paste0(
    data_fungi_mini@sam_data$Height,
    "_", data_fungi_mini@sam_data$Time
  )
  merge_samples2(data_fungi_mini, "H_T")
  hill_curves_pq(data_fungi_mini, "H_T", color_fac = "Time", nperm = 9)
}
#> 17 were discarded due to NA in variables present in formula.
#> At least one sample name start with a zero.
#>     That can be a problem for some phyloseq functions such as
#>     plot_bar and psmelt.
#> 47 were discarded due to NA in variables present in formula.
#> 47 were discarded due to NA in variables present in formula.
#> 'nperm' >= set of all permutations: complete enumeration.
#> Set of permutations < 'minperm'. Generating entire set.
#> 47 were discarded due to NA in variables present in formula.
#> 'nperm' >= set of all permutations: complete enumeration.
#> Set of permutations < 'minperm'. Generating entire set.
```
