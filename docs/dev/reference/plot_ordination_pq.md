# A wrapper of plot_ordination with vegan distance matrix

A wrapper of plot_ordination with vegan distance matrix

## Usage

``` r
plot_ordination_pq(
  physeq,
  method = "robust.aitchison",
  ordination_method = "NMDS",
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

- ...:

  Additional arguments passed on to phyloseq::plot_ordination()

## Value

A ggplot2 object

## Details

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wrapper of
[`phyloseq::plot_ordination()`](https://rdrr.io/pkg/phyloseq/man/plot_ordination.html)
to use aitchison and robust.aitchison distances from vegan package.

## Author

Adrien Taudière

## Examples

``` r
library(patchwork)
plot_ordination_pq(data_fungi_mini, method = "robust.aitchison", color = "Height") +
  plot_ordination_pq(data_fungi_mini, method = "bray", color = "Height")
#> Taxa are now in columns.
#> Run 0 stress 0.115973 
#> Run 1 stress 0.1183097 
#> Run 2 stress 0.1331099 
#> Run 3 stress 0.1242026 
#> Run 4 stress 0.1329001 
#> Run 5 stress 0.1298104 
#> Run 6 stress 0.1375308 
#> Run 7 stress 0.1381297 
#> Run 8 stress 0.1247712 
#> Run 9 stress 0.1183481 
#> Run 10 stress 0.1182422 
#> Run 11 stress 0.143682 
#> Run 12 stress 0.1247717 
#> Run 13 stress 0.1323693 
#> Run 14 stress 0.1183481 
#> Run 15 stress 0.1159623 
#> ... New best solution
#> ... Procrustes: rmse 0.001484625  max resid 0.01705688 
#> Run 16 stress 0.1266064 
#> Run 17 stress 0.1251491 
#> Run 18 stress 0.1400862 
#> Run 19 stress 0.1300025 
#> Run 20 stress 0.1782625 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     16: stress ratio > sratmax
#>      4: scale factor of the gradient < sfgrmin
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the phyloseq package.
#>   Please report the issue at <https://github.com/joey711/phyloseq/issues>.
#> Taxa are now in columns.
#> Run 0 stress 0.1894383 
#> Run 1 stress 0.1874393 
#> ... New best solution
#> ... Procrustes: rmse 0.07074037  max resid 0.256976 
#> Run 2 stress 0.185181 
#> ... New best solution
#> ... Procrustes: rmse 0.05747006  max resid 0.3051838 
#> Run 3 stress 0.1902199 
#> Run 4 stress 0.1879452 
#> Run 5 stress 0.191006 
#> Run 6 stress 0.1910298 
#> Run 7 stress 0.1854299 
#> ... Procrustes: rmse 0.02824792  max resid 0.1207013 
#> Run 8 stress 0.1870153 
#> Run 9 stress 0.1869737 
#> Run 10 stress 0.1870639 
#> Run 11 stress 0.1874451 
#> Run 12 stress 0.1896332 
#> Run 13 stress 0.1876277 
#> Run 14 stress 0.1897254 
#> Run 15 stress 0.1897057 
#> Run 16 stress 0.1911094 
#> Run 17 stress 0.1876672 
#> Run 18 stress 0.1896473 
#> Run 19 stress 0.1925673 
#> Run 20 stress 0.1850743 
#> ... New best solution
#> ... Procrustes: rmse 0.02856534  max resid 0.1295053 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     15: no. of iterations >= maxit
#>      5: stress ratio > sratmax
```
