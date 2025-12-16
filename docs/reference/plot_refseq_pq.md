# Plot the nucleotide proportion of references sequences

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

It is a wrapper of the function
[`plot_refseq_extremity_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_refseq_extremity_pq.md).
See ?plot_refseq_extremity_pq for more examples.

If `hill_scale` is not null, Hill diversity number are used to represent
the distribution of the diversity (equitability) along the sequences.

## Usage

``` r
plot_refseq_pq(
  physeq,
  hill_scales = NULL,
  first_n = min(Biostrings::width(physeq@refseq)),
  last_n = NULL,
  min_width = first_n
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- hill_scales:

  (vector) A vector defining the Hill number wanted. Set to NULL if you
  don't want to plot Hill diversity metrics.

- first_n:

  (int, default 10) The number of nucleotides to plot the 5' extremity.

- last_n:

  (int, default 10) The number of nucleotides to plot the 3' extremity.

- min_width:

  (int, default 0) Select only the sequences from physeq@refseq with
  using a minimum length threshold. If `first_n` is superior to the
  minimum length of the references sequences, you must use min_width to
  filter out the narrower sequences

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

plot_refseq_pq(data_fungi, hill_scales = c(2), first_n = 300)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 483
#> Number of kept ASV 937
#> Number of kept samples 185
#> Warning: Removed 492 rows containing missing values or values outside the scale range
#> (`geom_point()`).
#> Warning: Removed 123 rows containing missing values or values outside the scale range
#> (`geom_line()`).

```
