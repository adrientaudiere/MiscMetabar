# Calculate hill number and compute Tuckey post-hoc test

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that, by default, this function use a sqrt of the read numbers in
the linear model in order to correct for uneven sampling depth.

## Usage

``` r
hill_tuckey_pq(
  physeq,
  modality,
  hill_scales = c(0, 1, 2),
  silent = TRUE,
  correction_for_sample_size = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- modality:

  (required) the variable to test

- hill_scales:

  (a vector of integer) The list of q values to compute the hill number
  H^q. If Null, no hill number are computed. Default value compute the
  Hill number 0 (Species richness), the Hill number 1 (exponential of
  Shannon Index) and the Hill number 2 (inverse of Simpson Index).

- silent:

  (logical) If TRUE, no message are printing.

- correction_for_sample_size:

  (logical, default TRUE) This function use a sqrt of the read numbers
  in the linear model in order to correct for uneven sampling depth.

## Value

A ggplot2 object

## Author

Adrien Taudière

## Examples

``` r
data("GlobalPatterns", package = "phyloseq")
GlobalPatterns@sam_data[, "Soil_logical"] <-
  ifelse(GlobalPatterns@sam_data[, "SampleType"] == "Soil", "Soil", "Not Soil")
hill_tuckey_pq(GlobalPatterns, "Soil_logical")
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’

hill_tuckey_pq(GlobalPatterns, "Soil_logical", hill_scales = 1:2)
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
```
