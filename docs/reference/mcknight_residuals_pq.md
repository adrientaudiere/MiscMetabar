# Depth-robust alpha diversity residuals (McKnight / Mikryukov)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Computes the residuals from a linear regression of
\\\log(\text{richness})\\ against \\\log(\text{sequencing depth})\\ as a
depth-robust alpha diversity metric (Mikryukov et al. 2023; McKnight et
al. 2018,
[doi:10.5061/dryad.tn8qs35](https://doi.org/10.5061/dryad.tn8qs35) ).
This avoids discarding data through rarefaction.

## Usage

``` r
mcknight_residuals_pq(physeq, add_to_sam_data = TRUE)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- add_to_sam_data:

  (logical, default `TRUE`) if `TRUE`, a column `mcknight_residuals` is
  added to `sample_data(physeq)` and the augmented phyloseq object is
  returned. If `FALSE`, the numeric residuals vector is returned.

## Value

Either a phyloseq object with an augmented `sample_data` (default) or a
named numeric vector of residuals.

## Author

Adrien Taudière

## Examples

``` r
data_f_res <- mcknight_residuals_pq(data_fungi_mini)
head(sample_data(data_f_res)$mcknight_residuals)
#> [1]  0.3273090 -0.1430749  0.5681529  0.7310108 -0.7698491  0.4374109
```
