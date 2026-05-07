# Discard legend in ggplot2

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A more memorable shortcut for theme(legend.position = "none").

## Usage

``` r
no_legend()
```

## Value

A ggplot2 object

## Author

Adrien Taudière

## Examples

``` r
plot_refseq_pq(data_fungi)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1
#> Number of kept ASV 1419
#> Number of kept samples 185
#> Warning: Removed 688 rows containing missing values or values outside the scale range
#> (`geom_point()`).

plot_refseq_pq(data_fungi) + no_legend()
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1
#> Number of kept ASV 1419
#> Number of kept samples 185
#> Warning: Removed 688 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```
