# Calculate hill number and compute Tuckey post-hoc test

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that, by default, this function use a sqrt of the read numbers in
the linear model in order to correct for uneven sampling depth.

## Usage

``` r
hill_tuckey_pq(
  physeq,
  modality,
  q = c(0, 1, 2),
  hill_scales = lifecycle::deprecated(),
  silent = TRUE,
  correction_for_sample_size = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- modality:

  (required) the variable to test

- q:

  (numeric vector) Hill diversity orders to compute (q values). Default
  computes Hill number 0 (species richness), Hill number 1 (exponential
  of Shannon index) and Hill number 2 (inverse of Simpson index).
  Formerly `hill_scales`. Hill numbers are more appropriate in DNA
  metabarcoding studies when `q > 0` (Alberdi & Gilbert, 2019;
  Calderón-Sanou et al., 2019).

- hill_scales:

  **\[deprecated\]** Use `q` instead.

- silent:

  (logical) If TRUE, no message are printing.

- correction_for_sample_size:

  (logical, default TRUE) This function use a sqrt of the read numbers
  in the linear model in order to correct for uneven sampling depth.

- ...:

  Additional arguments passed to
  [`divent_hill_matrix_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/divent_hill_matrix_pq.md)
  and hence to
  [`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html)
  (e.g. `estimator = "naive"` to match vegan-style results).

## Value

A ggplot2 object

## References

Alberdi, A., & Gilbert, M. T. P. (2019). A guide to the application of
Hill numbers to DNA-based diversity analyses. *Molecular Ecology
Resources*.
[doi:10.1111/1755-0998.13014](https://doi.org/10.1111/1755-0998.13014)

Calderón-Sanou, I., Münkemüller, T., Boyer, F., Zinger, L., & Thuiller,
W. (2019). From environmental DNA sequences to ecological conclusions:
How strong is the influence of methodological choices? *Journal of
Biogeography*, 47.
[doi:10.1111/jbi.13681](https://doi.org/10.1111/jbi.13681)

## Author

Adrien Taudière

## Examples

``` r
data("GlobalPatterns", package = "phyloseq")
GlobalPatterns@sam_data[, "Soil_logical"] <-
  ifelse(GlobalPatterns@sam_data[, "SampleType"] == "Soil", "Soil", "Not Soil")
hill_tuckey_pq(GlobalPatterns, "Soil_logical")

hill_tuckey_pq(GlobalPatterns, "Soil_logical", q = 1:2)
```
