# Plot taxonomic distribution across 3 taxonomic levels and optionally one sample factor

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that lvl3 need to be nested in lvl2 which need to be nested in lvl1

## Usage

``` r
multitax_bar_pq(
  physeq,
  lvl1,
  lvl2,
  lvl3,
  fact = NULL,
  nb_seq = TRUE,
  log10trans = TRUE
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

  (required) Name of the second (middle) taxonomic rank of interest

- lvl3:

  (required) Name of the first (lower) taxonomic rank of interest

- fact:

  Name of the factor to cluster samples by modalities. Need to be in
  `physeq@sam_data`. If not set, the taxonomic distribution is plot for
  all samples together.

- nb_seq:

  (logical; default TRUE) If set to FALSE, only the number of ASV is
  count. Concretely, physeq otu_table is transformed in a binary
  otu_table (each value different from zero is set to one)

- log10trans:

  (logical, default TRUE) If TRUE, the number of sequences (or ASV if
  nb_seq = FALSE) is log10 transformed.

## Value

A ggplot2 object

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("ggh4x")) {
  multitax_bar_pq(data_fungi_sp_known, "Phylum", "Class", "Order", "Time")
  multitax_bar_pq(data_fungi_sp_known, "Phylum", "Class", "Order")
  multitax_bar_pq(data_fungi_sp_known, "Phylum", "Class", "Order",
    nb_seq = FALSE, log10trans = FALSE
  )
}
#> Loading required namespace: ggh4x

# }
```
