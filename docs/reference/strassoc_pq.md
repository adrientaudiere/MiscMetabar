# Strength of species-group associations for a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`indicspecies::strassoc()`](https://emf-creaf.github.io/indicspecies/reference/strassoc.html)
function in the case of a `physeq` object. It computes the strength of
the association between each taxon and the groups defined by `fact`,
using one of the association indices provided by the `indicspecies`
package (e.g. `IndVal`, the phi coefficient `r`, the specificity
component `A`, or the fidelity component `B`). It complements
[`multipatt_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multipatt_pq.md)
which tests the statistical significance of those associations.

## Usage

``` r
strassoc_pq(
  physeq,
  fact,
  func = "IndVal.g",
  nboot_ci = NULL,
  alpha_ci = 0.05,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor in `physeq@sam_data` defining the groups
  of samples.

- func:

  (chr, default "IndVal.g") The association index to compute. See
  `?indicspecies::strassoc()` for the full list. Common values are:

  - `"IndVal"` / `"IndVal.g"`: Dufrêne-Legendre indicator value (the
    `.g` variant equalizes group sizes),

  - `"r"` / `"r.g"`: the point-biserial / phi coefficient of
    association,

  - `"A"`: specificity (positive predictive value) component,

  - `"B"`: fidelity (sensitivity) component.

- nboot_ci:

  (int, default NULL) Number of bootstrap replicates used to estimate
  confidence intervals. When NULL (default), no bootstrap is performed
  and a single estimate is returned per taxon and group.

- alpha_ci:

  (float, default 0.05) Error level for the bootstrap confidence
  intervals (only used when `nboot_ci` is not NULL).

- ...:

  Additional arguments passed on to
  [`indicspecies::strassoc()`](https://emf-creaf.github.io/indicspecies/reference/strassoc.html).

## Value

- When `nboot_ci` is NULL: a `tibble` with one row per taxon, a `taxon`
  column and one column per group containing the association value.

- When `nboot_ci` is set: the list returned by
  [`indicspecies::strassoc()`](https://emf-creaf.github.io/indicspecies/reference/strassoc.html)
  with the `lowerCI`, `stat` and `upperCI` matrices.

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`indicspecies::strassoc()`](https://emf-creaf.github.io/indicspecies/reference/strassoc.html)
if you use this function.

## See also

[`multipatt_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multipatt_pq.md),
[`indicspecies::strassoc()`](https://emf-creaf.github.io/indicspecies/reference/strassoc.html)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("indicspecies")) {
  strassoc_pq(
    subset_samples(data_fungi_mini, !is.na(Height)),
    fact = "Height",
    func = "IndVal.g"
  )
}
#> Taxa are now in columns.
#> # A tibble: 45 × 4
#>    taxon   High    Low  Middle
#>    <chr>  <dbl>  <dbl>   <dbl>
#>  1 ASV7  0.519  0.157  0.268  
#>  2 ASV8  0.553  0.239  0.474  
#>  3 ASV12 0.178  0.391  0.157  
#>  4 ASV18 0.508  0.0984 0.0939 
#>  5 ASV25 0.459  0.0229 0.0218 
#>  6 ASV26 0.0173 0.105  0.304  
#>  7 ASV27 0.0576 0.168  0.00346
#>  8 ASV29 0.0340 0.229  0.266  
#>  9 ASV32 0.266  0.0169 0      
#> 10 ASV34 0.154  0.266  0.387  
#> # ℹ 35 more rows
# }
if (FALSE) { # \dontrun{
if (requireNamespace("indicspecies")) {
  strassoc_pq(
    subset_samples(data_fungi_mini, !is.na(Height)),
    fact = "Height",
    func = "A"
  )
  strassoc_pq(
    subset_samples(data_fungi_mini, !is.na(Height)),
    fact = "Height",
    func = "IndVal.g",
    nboot_ci = 99
  )
}
} # }
```
