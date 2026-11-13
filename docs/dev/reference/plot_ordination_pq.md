# A wrapper of plot_ordination with vegan distance matrix

A wrapper of plot_ordination with vegan distance matrix

## Usage

``` r
plot_ordination_pq(
  physeq,
  method = "robust.aitchison",
  ordination_method = "NMDS",
  add_envfit = FALSE,
  envfit_fact = NULL,
  envfit_pval = 1,
  envfit_arrow_mult = 1,
  envfit_text_size = 3,
  envfit_col = "black",
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- method:

  (string, default "robust.aitchison") The distance method to use from
  vegan::vegdist(). See ?vegan::vegdist for more details.

- ordination_method:

  (string, default "NMDS") The ordination method to use in
  phyloseq::ordinate(). See ?phyloseq::ordinate for more details.

- add_envfit:

  (logical, default FALSE) If TRUE, overlay environmental vectors (and
  factor centroids) fitted to the ordination using
  [`vegan::envfit()`](https://vegandevs.github.io/vegan/reference/envfit.html).
  Continuous variables are shown as arrows and factor variables as
  labelled centroids.

- envfit_fact:

  (default NULL) A character vector of variable names from
  `physeq@sam_data` to fit onto the ordination. When NULL (default) and
  `add_envfit = TRUE`, every variable of `physeq@sam_data` is used;
  passing an explicit subset is recommended, as fitting identifier or
  high-cardinality character columns produces a cluttered overlay.

- envfit_pval:

  (float, default 1) Only environmental variables with an
  [`vegan::envfit()`](https://vegandevs.github.io/vegan/reference/envfit.html)
  p-value lower than or equal to this threshold are drawn. The
  default (1) draws all of them.

- envfit_arrow_mult:

  (float, default 1) A multiplier applied to the environmental arrow
  lengths, to scale them for readability.

- envfit_text_size:

  (float, default 3) Text size for the environmental variable labels.

- envfit_col:

  (default "black") Color of the environmental arrows, centroids and
  labels.

- ...:

  Additional arguments passed on to phyloseq::plot_ordination()

## Value

A ggplot2 object

## Details

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wrapper of
[`phyloseq::plot_ordination()`](https://rdrr.io/pkg/phyloseq/man/plot_ordination.html)
to use aitchison and robust.aitchison distances from vegan package.

## See also

[`vegan::envfit()`](https://vegandevs.github.io/vegan/reference/envfit.html),
[`phyloseq::plot_ordination()`](https://rdrr.io/pkg/phyloseq/man/plot_ordination.html)

## Author

Adrien Taudière

## Examples

``` r
library(patchwork)
plot_ordination_pq(data_fungi_mini, method = "robust.aitchison", color = "Height") +
  plot_ordination_pq(data_fungi_mini, method = "bray", color = "Height")
#> Taxa are now in columns.
#> Run 0 stress 0.115973 
#> Run 1 stress 0.1330388 
#> Run 2 stress 0.1735524 
#> Run 3 stress 0.1430453 
#> Run 4 stress 0.1183481 
#> Run 5 stress 0.1389906 
#> Run 6 stress 0.1321663 
#> Run 7 stress 0.1795299 
#> Run 8 stress 0.1297006 
#> Run 9 stress 0.11883 
#> Run 10 stress 0.1266062 
#> Run 11 stress 0.1242026 
#> Run 12 stress 0.1204173 
#> Run 13 stress 0.1187081 
#> Run 14 stress 0.1561721 
#> Run 15 stress 0.1339447 
#> Run 16 stress 0.1183481 
#> Run 17 stress 0.142358 
#> Run 18 stress 0.1183481 
#> Run 19 stress 0.1554865 
#> Run 20 stress 0.1159623 
#> ... New best solution
#> ... Procrustes: rmse 0.001485021  max resid 0.01705947 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>      1: no. of iterations >= maxit
#>     10: stress ratio > sratmax
#>      9: scale factor of the gradient < sfgrmin
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the phyloseq package.
#>   Please report the issue at <https://github.com/joey711/phyloseq/issues>.
#> Taxa are now in columns.
#> Run 0 stress 0.1894178 
#> Run 1 stress 0.192511 
#> Run 2 stress 0.1856177 
#> ... New best solution
#> ... Procrustes: rmse 0.06395323  max resid 0.2137433 
#> Run 3 stress 0.1876421 
#> Run 4 stress 0.1889902 
#> Run 5 stress 0.1880454 
#> Run 6 stress 0.1870819 
#> Run 7 stress 0.191065 
#> Run 8 stress 0.1871421 
#> Run 9 stress 0.1914317 
#> Run 10 stress 0.1917667 
#> Run 11 stress 0.1925024 
#> Run 12 stress 0.187398 
#> Run 13 stress 0.1903712 
#> Run 14 stress 0.1873517 
#> Run 15 stress 0.1884903 
#> Run 16 stress 0.1888557 
#> Run 17 stress 0.1872554 
#> Run 18 stress 0.1905621 
#> Run 19 stress 0.189407 
#> Run 20 stress 0.1869003 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     16: no. of iterations >= maxit
#>      4: stress ratio > sratmax
```
