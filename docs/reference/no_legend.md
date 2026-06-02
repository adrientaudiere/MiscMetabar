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
plot_refseq_pq(data_fungi_mini)
#> Cleaning suppress 0 taxa (  ) and 1 sample(s) ( W26-001-B_S165_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 1
#> Number of kept ASV 44
#> Number of kept samples 136
#> Warning: Removed 336 rows containing missing values or values outside the scale range
#> (`geom_point()`).

plot_refseq_pq(data_fungi_mini) + no_legend()
#> Cleaning suppress 0 taxa (  ) and 1 sample(s) ( W26-001-B_S165_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 1
#> Number of kept ASV 44
#> Number of kept samples 136
#> Warning: Removed 336 rows containing missing values or values outside the scale range
#> (`geom_point()`).
```
