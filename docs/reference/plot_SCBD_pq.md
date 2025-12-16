# Plot species contributions to beta diversity (SCBD) of samples

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
function in the case of `physeq` object.

## Usage

``` r
plot_SCBD_pq(
  physeq,
  tax_level = "Taxa",
  tax_col = "Order",
  min_SCBD = 0.01,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- tax_level:

  Taxonomic level to used in y axis

- tax_col:

  Taxonomic level to colored points

- min_SCBD:

  (default 0.01) the minimum SCBD value to plot the taxa

- ...:

  Additional arguments passed on to
  [`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
  function

## Value

A ggplot2 object build with the package patchwork

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to `vegan::beta.div()` if you use this function.

## See also

[LCBD_pq](https://adrientaudiere.github.io/MiscMetabar/reference/LCBD_pq.md),
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)

## Author

Adrien Taudière

## Examples

``` r
data(data_fungi)
if (requireNamespace("adespatial")) {
  plot_SCBD_pq(data_fungi) +
    geom_text(aes(label = paste(Genus, Species)), hjust = 1, vjust = 2) +
    xlim(c(0, NA))
}

# \donttest{
if (requireNamespace("adespatial")) {
  plot_SCBD_pq(data_fungi, tax_level = "Class", tax_col = "Phylum", min_SCBD = 0) +
    geom_jitter()
}

# }
```
