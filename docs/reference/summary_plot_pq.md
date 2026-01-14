# Summarize a [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object using a plot.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of a phyloseq object.

## Usage

``` r
summary_plot_pq(
  physeq,
  add_info = TRUE,
  min_seq_samples = 500,
  clean_pq = TRUE,
  text_size = 1,
  text_size_info = 1
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- add_info:

  Does the bottom down corner contain extra informations?

- min_seq_samples:

  (int): Used only when add_info is set to true to print the number of
  samples with less sequences than this number.

- clean_pq:

  (logical): Does the phyloseq object is cleaned using the
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
  function?

- text_size:

  (Num, default 1) A size factor to expand or minimize text size.

- text_size_info:

  (Num, default 1) A size factor to expand or minimize text size for
  extra informations.

## Value

A ggplot2 object

## Examples

``` r
summary_plot_pq(data_fungi)

summary_plot_pq(data_fungi, add_info = FALSE) + scale_fill_viridis_d()
#> Scale for fill is already present.
#> Adding another scale for fill, which will replace the existing scale.

if (requireNamespace("patchwork")) {
  (summary_plot_pq(data_fungi, text_size = 0.5, text_size_info = 0.6) +
    summary_plot_pq(data_fungi_mini, text_size = 0.5, text_size_info = 0.6)) /
    (summary_plot_pq(data_fungi_sp_known, text_size = 0.5, text_size_info = 0.6) +
      summary_plot_pq(subset_taxa(data_fungi_sp_known, Phylum == "Ascomycota"),
        text_size = 0.5, text_size_info = 0.6
      ))
}
#> Cleaning suppress 0 taxa and 1 samples.
#> Cleaning suppress 0 taxa and 2 samples.
```
