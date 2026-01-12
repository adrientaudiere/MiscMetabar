# Run ANCOMBC2 on phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html)
function

## Usage

``` r
ancombc_pq(physeq, fact, levels_fact = NULL, tax_level = "Class", ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor in `physeq@sam_data` used to plot
  different lines

- levels_fact:

  (default NULL) The order of the level in the factor. Used for reorder
  levels and select levels (filter out levels not present en
  levels_fact)

- tax_level:

  The taxonomic level passed on to
  [`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html)

- ...:

  Additional arguments passed on to
  [`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html)
  function.

## Value

The result of
[`ANCOMBC::ancombc2()`](https://rdrr.io/pkg/ANCOMBC/man/ancombc2.html)
function

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
  res_height <- ancombc_pq(
    data_fungi_mini,
    fact = "Height",
    levels_fact = c("Low", "High"),
    verbose = TRUE
  )

  ggplot(
    res_height$res,
    aes(
      y = reorder(taxon, lfc_HeightHigh),
      x = lfc_HeightHigh,
      color = diff_HeightHigh
    )
  ) +
    geom_vline(xintercept = 0) +
    geom_segment(aes(
      xend = 0, y = reorder(taxon, lfc_HeightHigh),
      yend = reorder(taxon, lfc_HeightHigh)
    ), color = "darkgrey") +
    geom_point()

  res_time <- ancombc_pq(
    data_fungi_mini,
    fact = "Time",
    levels_fact = c("0", "15"),
    tax_level = "Family",
    verbose = TRUE
  )
}
} # }
```
