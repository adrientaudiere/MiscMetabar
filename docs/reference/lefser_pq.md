# Run LEfSe on a phyloseq object

Run LEfSe on a phyloseq object

## Usage

``` r
lefser_pq(
  physeq,
  bifactor = NULL,
  modalities = NULL,
  compute_relativeAb = TRUE,
  by_clade = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- bifactor:

  (required) The name of a column present in the `@sam_data` slot of the
  physeq object. Must be a character vector or a factor.

- modalities:

  (default NULL) A vector of modalities to keep in the analysis. If
  NULL, all modalities present in bifactor are kept. Note that only two
  modalities are allowed. @param compute_relativeAb (logical, default
  TRUE) Do we compute relative abundance before running LEfSe? @param
  by_clade (logical, default FALSE) Do we use the lefserClades function
  (which test for different depth in the taxonomic classification) or
  the lefser function (taxa-level)? @param ... Additional arguments
  passed on to
  [`lefser::lefser()`](https://rdrr.io/pkg/lefser/man/lefser.html)

## Value

The result of lefser::lefser() or lefser::lefserClades()

## Details

It is a wrapper of the
[`lefser::lefser()`](https://rdrr.io/pkg/lefser/man/lefser.html) and
[`lefser::lefserClades()`](https://rdrr.io/pkg/lefser/man/lefserClades.html)
functions.
[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Author

Adrien Taudière

## Examples

``` r
res_lefse <- lefser_pq(data_fungi,
  bifactor = "Height",
  modalities = c("Low", "High")
)
#> Error in relativeAb(physeq_ts): could not find function "relativeAb"
lefser::lefserPlot(res_lefse)
#> Error: object 'res_lefse' not found
```
