# Visualization of a collection of couples of samples for comparison

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This allow to plot all the possible
[`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)
combination using one factor.

## Usage

``` r
multi_biplot_pq(physeq, split_by = NULL, pairs = NULL, na_remove = TRUE, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- split_by:

  (required if pairs is NULL) the name of the factor to make all
  combination of couples of values

- pairs:

  (required if split_by is NULL) the name of the factor in
  physeq@sam_data\` slot to make plot by pairs of samples. Each level
  must be present only two times. Note that if you set pairs, you also
  must set fact arguments to passed on to
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md).

- na_remove:

  (logical, default TRUE) if TRUE remove all the samples with NA in the
  `split_by` variable of the `physeq@sam_data` slot

- ...:

  Other parameters passed on to
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)

## Value

a list of ggplot object

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
data_fungi_abun <- subset_taxa_pq(data_fungi_mini,
  taxa_sums(data_fungi_mini) > 1000)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 0
#> Number of kept ASV 45
#> Number of kept samples 137
p <- multi_biplot_pq(data_fungi_abun, "Height")
#> 47 were discarded due to NA in variables present in formula.
#> Cleaning suppress 5 taxa and 0 samples.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
#> Cleaning suppress 4 taxa and 0 samples.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
#> Cleaning suppress 7 taxa and 0 samples.
#> Scale for y is already present.
#> Adding another scale for y, which will replace the existing scale.
lapply(p, print)



#> $`High-Low`

#> 
#> $`High-Middle`

#> 
#> $`Low-Middle`

#> 
# }
```
