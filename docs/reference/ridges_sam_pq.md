# Ridges plot of sample distribution across taxa

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of distribution of samples across taxa using
ridges. This is the sample-centric counterpart of
[`ridges_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ridges_pq.md):
each ridge represents a taxon (at `tax_level`) and the x-axis shows the
abundance distribution across samples, optionally colored by a sample
factor.

## Usage

``` r
ridges_sam_pq(
  physeq,
  fact,
  nb_seq = TRUE,
  log10trans = TRUE,
  tax_level = "Class",
  type = "density",
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor in `physeq@sam_data` used to color the
  ridges

- nb_seq:

  (logical; default TRUE) If set to FALSE, only the number of samples is
  counted. Concretely, physeq `otu_table` is transformed in a binary
  `otu_table` (each value different from zero is set to one)

- log10trans:

  (logical, default TRUE) If TRUE, the abundance is log10 transformed.

- tax_level:

  The taxonomic level used for grouping taxa on the y-axis

- type:

  Either "density" (the default) or "ecdf" to plot a cumulative version
  using
  [`ggplot2::stat_ecdf()`](https://ggplot2.tidyverse.org/reference/stat_ecdf.html)

- ...:

  Other params passed on to
  [`ggridges::geom_density_ridges()`](https://wilkelab.org/ggridges/reference/geom_density_ridges.html)

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot
with ridges representing the distribution of samples for each taxon

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("ggridges")) {
  ridges_sam_pq(data_fungi_mini, "Height", alpha = 0.5,
    log10trans = FALSE, tax_level = "Genus") +
  xlim(c(0, 1000))
}
#> Loading required namespace: ggridges
#> Scale for x is already present.
#> Adding another scale for x, which will replace the existing scale.
#> Picking joint bandwidth of 36.4
#> Warning: Removed 99 rows containing non-finite outside the scale range
#> (`stat_density_ridges()`).

# \donttest{
if (requireNamespace("ggridges")) {
  ridges_sam_pq(data_fungi_mini, "Height", alpha = 0.5, scale = 0.9)
  ridges_sam_pq(data_fungi_mini, "Height",
    alpha = 0.5, scale = 0.9,
    type = "ecdf"
  )
}

# }
```
