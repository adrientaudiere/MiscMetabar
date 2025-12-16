# Plot information about Guild from tax_table slot previously created with [`add_funguild_info()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_funguild_info.md)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical function.

## Usage

``` r
plot_guild_pq(physeq, levels_order = NULL, clean_pq = TRUE, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- levels_order:

  (Default NULL) A character vector to reorder the levels of guild. See
  examples.

- clean_pq:

  (logical, default TRUE): Does the phyloseq object is cleaned using the
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
  function?

- ...:

  Other params for be passed on to
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
  function

## Value

A ggplot2 object

## See also

[`add_funguild_info()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_funguild_info.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# to avoid bug in CRAN when internet is not available
if (requireNamespace("httr")) {
  d_fung_mini <- add_funguild_info(data_fungi_mini,
    taxLevels = c(
      "Domain",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )
  )
  sort(table(d_fung_mini@tax_table[, "guild"]), decreasing = TRUE)

  p <- plot_guild_pq(d_fung_mini)
  if (requireNamespace("patchwork")) {
    (plot_guild_pq(subset_samples(d_fung_mini, Height == "Low"),
      levels_order = p$data$Guild[order(p$data$nb_seq)]
    ) + theme(legend.position = "none")) +
      (plot_guild_pq(subset_samples(d_fung_mini, Height == "High"),
        levels_order = p$data$Guild[order(p$data$nb_seq)]
      ) + ylab("") + theme(axis.text.y = element_blank()))
  }
}
} # }
```
