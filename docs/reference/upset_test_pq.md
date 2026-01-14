# Test for differences between intersections

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

See
[`upset_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/upset_pq.md)
to plot upset. There is a bug with ggplot2 \>= 4.0.0. See issue
<https://github.com/krassowski/complex-upset/issues/213> for more
details.

## Usage

``` r
upset_test_pq(
  physeq,
  fact,
  var_to_test = "OTU",
  min_nb_seq = 0,
  na_remove = TRUE,
  numeric_fonction = sum,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor to cluster samples by modalities. Need
  to be in `physeq@sam_data`.

- var_to_test:

  (default c("OTU")) : a vector of column present in the tax_table slot
  from the physeq object

- min_nb_seq:

  minimum number of sequences by OTUs by samples to take into count this
  OTUs in this sample. For example, if min_nb_seq=2,each value of 2 or
  less in the OTU table will not count in the venn diagram

- na_remove:

  : if TRUE (the default), NA values in fact are removed if FALSE, NA
  values are set to "NA"

- numeric_fonction:

  (default : sum) the function for numeric vector useful only for
  complex plot (see examples)

- ...:

  Additional arguments passed on to the
  [`ComplexUpset::upset_test()`](https://krassowski.github.io/complex-upset/reference/upset_test.html)

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot

## See also

[`upset_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/upset_pq.md)

## Author

Adrien Taudière

## Examples

``` r
data(data_fungi)
if (requireNamespace("ComplexUpset")) {
  upset_test_pq(data_fungi, "Height", var_to_test = c("OTU", "Class", "Guild"))
  upset_test_pq(data_fungi, "Time")
}
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 0
#> Number of kept ASV 1420
#> Number of kept samples 185
#> [1] "OTU, Class, Guild differ significantly between intersections"
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 0
#> Number of kept ASV 1420
#> Number of kept samples 185
#> [1] "OTU differ significantly between intersections"
#>     variable     p.value statistic                         test         fdr
#> OTU      OTU 0.002217898   33.7837 Kruskal-Wallis rank sum test 0.002217898
```
