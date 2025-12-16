# Plot ANCOMBC2 result for phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of ANCOMBC2 result.

## Usage

``` r
plot_ancombc_pq(
  physeq,
  ancombc_res,
  filter_passed = TRUE,
  filter_diff = TRUE,
  min_abs_lfc = 0,
  tax_col = "Genus",
  tax_label = "Species",
  add_marginal_vioplot = TRUE,
  add_label = TRUE,
  add_hline_cut_lfc = NULL
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

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

  (integer, default 0) Minimum absolute value to filter results based on
  Log Fold Change. For ex. a value of 1 filter out taxa for which the
  abundance in a given level of the modalty is not at least the double
  of the abundance in the other level.

- tax_col:

  The taxonomic level (must be present in `tax_table` slot) to color the
  points

- tax_label:

  The taxonomic level (must be present in `tax_table` slot) to add label

- add_marginal_vioplot:

  (logical, default TRUE) Do we add a marginal vioplot representing all
  the taxa lfc from ancombc_res.

- add_label:

  (logical, default TRUE) Do we add a label?

- add_hline_cut_lfc:

  (logical, default NULL) Do we add two horizontal lines when
  min_abs_lfc is set (different from zero)?

## Value

A ggplot2 object. If add_marginal_vioplot is TRUE, this is a patchworks
of plot made using
[`patchwork::plot_layout()`](https://patchwork.data-imaginist.com/reference/plot_layout.html).

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html)
if you use this function.

## Author

Adrien Taudière

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

  plot_ancombc_pq(data_fungi_mini, res_time,
    filter_passed = FALSE,
    tax_label = "Genus", tax_col = "Order"
  )
  plot_ancombc_pq(data_fungi_mini, res_time, tax_col = "Genus")
  plot_ancombc_pq(data_fungi_mini, res_time,
    filter_passed = FALSE,
    filter_diff = FALSE, tax_col = "Family", add_label = FALSE
  )
}
} # }
```
