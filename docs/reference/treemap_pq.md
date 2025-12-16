# Plot treemap of 2 taxonomic levels

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that lvl2need to be nested in lvl1

## Usage

``` r
treemap_pq(
  physeq,
  lvl1,
  lvl2,
  nb_seq = TRUE,
  log10trans = TRUE,
  plot_legend = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- lvl1:

  (required) Name of the first (higher) taxonomic rank of interest

- lvl2:

  (required) Name of the second (lower) taxonomic rank of interest

- nb_seq:

  (logical; default TRUE) If set to FALSE, only the number of ASV is
  count. Concretely, physeq otu_table is transformed in a binary
  otu_table (each value different from zero is set to one)

- log10trans:

  (logical, default TRUE) If TRUE, the number of sequences (or ASV if
  nb_seq = FALSE) is log10 transformed.

- plot_legend:

  (logical, default FALSE) If TRUE, plot che legend of color for lvl 1

- ...:

  Additional arguments passed on to
  [`treemapify::geom_treemap()`](https://wilkox.org/treemapify/reference/geom_treemap.html)
  function.

## Value

A ggplot2 object

## Author

Adrien Taudière

## Examples

``` r
data(data_fungi_sp_known)
if (requireNamespace("treemapify")) {
  treemap_pq(
    clean_pq(subset_taxa(
      data_fungi_sp_known,
      Phylum == "Basidiomycota"
    )),
    "Order", "Class",
    plot_legend = TRUE
  )
}
#> Loading required namespace: treemapify
#> Cleaning suppress 0 taxa and 5 samples.

# \donttest{
if (requireNamespace("treemapify")) {
  treemap_pq(
    clean_pq(subset_taxa(
      data_fungi_sp_known,
      Phylum == "Basidiomycota"
    )),
    "Order", "Class",
    log10trans = FALSE
  )
  treemap_pq(
    clean_pq(subset_taxa(
      data_fungi_sp_known,
      Phylum == "Basidiomycota"
    )),
    "Order", "Class",
    nb_seq = FALSE, log10trans = FALSE
  )
}
#> Cleaning suppress 0 taxa and 5 samples.
#> Cleaning suppress 0 taxa and 5 samples.

# }
```
