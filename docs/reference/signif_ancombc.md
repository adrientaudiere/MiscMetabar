# Filter ancombc_pq results

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Internally used in
[`plot_ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_ancombc_pq.md).

## Usage

``` r
signif_ancombc(
  ancombc_res,
  filter_passed = TRUE,
  filter_diff = TRUE,
  min_abs_lfc = 0
)
```

## Arguments

- ancombc_res:

  (required) the result of the ancombc_pq function For the moment only
  bimodal factors are possible.

- filter_passed:

  (logical, default TRUE) Do we filter using the column passed_ss? The
  passed_ss value is TRUE if the taxon passed the sensitivity analysis,
  i.e., adding different pseudo-counts to 0s would not change the
  results.

- filter_diff:

  (logical, default TRUE) Do we filter using the column diff? The diff
  value is TRUE if the taxon is significant (has q less than alpha)

- min_abs_lfc:

  (integer, default0) Minimum absolute value to filter results based on
  Log Fold Change. For ex. a value of 1 filter out taxa for which the
  abundance in a given level of the modality is not at least the double
  of the abundance in the other level.

## Value

A data.frame with the same number of columns than the `ancombc_res`
param but with less (or equal) numbers of rows

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html)
if you use this function.

## See also

[`ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ancombc_pq.md),
[`plot_ancombc_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_ancombc_pq.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("mia")) {
  data_fungi_mini@tax_table <- phyloseq::tax_table(cbind(
    data_fungi_mini@tax_table,
    "taxon" = taxa_names(data_fungi_mini)
  ))

  res_time <- ancombc_pq(
    data_fungi_mini,
    fact = "Time",
    levels_fact = c("0", "15"),
    tax_level = "taxon",
    verbose = TRUE
  )

  signif_ancombc(res_time)
}
} # }
```
