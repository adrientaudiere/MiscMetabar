# Plot and test local contributions to beta diversity (LCBD) of samples

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
function in the case of `physeq` object.

## Usage

``` r
plot_LCBD_pq(
  physeq,
  p_adjust_method = "BH",
  pval = 0.05,
  sam_variables = NULL,
  only_plot_significant = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- p_adjust_method:

  (chr, default "BH"): the method used to adjust p-value

- pval:

  (int, default 0.05): the value to determine the significance of LCBD

- sam_variables:

  A vector of variable names present in the `sam_data` slot to plot
  alongside the LCBD value

- only_plot_significant:

  (logical, default TRUE) Do we plot all LCBD values or only the
  significant ones

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
  plot_LCBD_pq(data_fungi_mini,
    nperm = 100, only_plot_significant = FALSE,
    pval = 0.2
  )
}
#> Taxa are now in columns.

# \donttest{
if (requireNamespace("adespatial")) {
  plot_LCBD_pq(data_fungi_mini,
    nperm = 100, only_plot_significant = TRUE,
    pval = 0.2
  )
  if (requireNamespace("patchwork")) {
    plot_LCBD_pq(data_fungi_mini,
      nperm = 100, only_plot_significant = FALSE,
      sam_variables = c("Time", "Height")
    )
    plot_LCBD_pq(data_fungi_mini,
      nperm = 100, only_plot_significant = TRUE, pval = 0.2,
      sam_variables = c("Time", "Height", "Tree_name")
    ) &
      theme(
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 10),
        axis.title.x = element_text(size = 6)
      )
  }
}
#> Taxa are now in columns.
#> Taxa are now in columns.
#> Taxa are now in columns.

# }
```
