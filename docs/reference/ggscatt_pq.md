# Scatterplot with marginal distributions and statistical results against Hill diversity of phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wrapper of function
[`ggstatsplot::ggscatterstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggscatterstats.html)
for object of class phyloseq and Hill number.

## Usage

``` r
ggscatt_pq(
  physeq,
  num_modality,
  hill_scales = c(0, 1, 2),
  rarefy_by_sample = FALSE,
  rngseed = FALSE,
  verbose = TRUE,
  one_plot = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- num_modality:

  (required) Name of the numeric column in `physeq@sam_data` to plot and
  test against hill number

- hill_scales:

  (a vector of integer) The list of q values to compute the hill number
  H^q. If Null, no hill number are computed. Default value compute the
  Hill number 0 (Species richness), the Hill number 1 (exponential of
  Shannon Index) and the Hill number 2 (inverse of Simpson Index).

- rarefy_by_sample:

  (logical, default FALSE) If TRUE, rarefy samples using
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  function.

- rngseed:

  (Optional). A single integer value passed to
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html),
  which is used to fix a seed for reproducibly random number generation
  (in this case, reproducibly random subsampling). If set to FALSE, then
  no fiddling with the RNG seed is performed, and it is up to the user
  to appropriately call set.seed beforehand to achieve reproducible
  results. Default is FALSE.

- verbose:

  (logical). If TRUE, print additional information.

- one_plot:

  (logical, default FALSE) If TRUE, return a unique plot with the three
  plot inside using the patchwork package.

- ...:

  Additional arguments passed on to
  [`ggstatsplot::ggscatterstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggscatterstats.html)
  function.

## Value

Either an unique ggplot2 object (if one_plot is TRUE) or a list of
ggplot2 plot for each hill_scales.

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`ggstatsplot::ggscatterstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggscatterstats.html)
if you use this function.

## See also

[`ggbetween_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggbetween_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("ggstatsplot")) {
  ggscatt_pq(data_fungi_mini, "Time", type = "non-parametric")
  ggscatt_pq(data_fungi_mini, "Time", hill_scales = 1:4, type = "parametric")
  ggscatt_pq(data_fungi_mini, "Sample_id",
    hill_scales = c(0, 0.5),
    one_plot = FALSE
  )
}
#> Taxa are now in columns.
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Taxa are now in columns.
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Taxa are now in columns.
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> [[1]]
#> `stat_xsidebin()` using `bins = 30`. Pick better value `binwidth`.
#> `stat_ysidebin()` using `bins = 30`. Pick better value `binwidth`.

#> 
#> [[2]]
#> `stat_xsidebin()` using `bins = 30`. Pick better value `binwidth`.
#> `stat_ysidebin()` using `bins = 30`. Pick better value `binwidth`.

#> 
```
