# Plot the distribution of sequences or ASV in one taxonomic levels

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of distribution of taxonomy, optionnaly across
a factor.

## Usage

``` r
tax_bar_pq(
  physeq,
  fact = "Sample",
  taxa = "Order",
  percent_bar = FALSE,
  nb_seq = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  Name of the factor to cluster samples by modalities. Need to be in
  `physeq@sam_data`.

- taxa:

  (default: 'Order') Name of the taxonomic rank of interest

- percent_bar:

  (default FALSE) If TRUE, the stacked bar fill all the space between 0
  and 1. It just set position = "fill" in the
  [`ggplot2::geom_bar()`](https://ggplot2.tidyverse.org/reference/geom_bar.html)
  function

- nb_seq:

  (logical; default TRUE) If set to FALSE, only the number of ASV is
  count. Concretely, physeq otu_table is transformed in a binary
  otu_table (each value different from zero is set to one)

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot
with bar representing the number of sequence en each taxonomic groups

## See also

[`plot_tax_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_tax_pq.md)
and
[`multitax_bar_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/multitax_bar_pq.md)

## Author

Adrien Taudière

## Examples

``` r
data_fungi_ab <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)
#> Cleaning suppress 0 taxa (  ) and 15 sample(s) ( BE9-006-B_S27_MERGED.fastq.gz / C21-NV1-M_S64_MERGED.fastq.gz / DJ2-008-B_S87_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / DY5-004-M_S98_MERGED.fastq.gz / E9-009-B_S100_MERGED.fastq.gz / E9-009-H_S101_MERGED.fastq.gz / N22-001-B_S129_MERGED.fastq.gz / O20-X-B_S139_MERGED.fastq.gz / O21-007-M_S144_MERGED.fastq.gz / R28-008-H_S159_MERGED.fastq.gz / R28-008-M_S160_MERGED.fastq.gz / W26-001-M_S167_MERGED.fastq.gz / Y29-007-H_S182_MERGED.fastq.gz / Y29-007-M_S183_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1385
#> Number of kept ASV 35
#> Number of kept samples 170
tax_bar_pq(data_fungi_ab) + theme(legend.position = "none")

# \donttest{
tax_bar_pq(data_fungi_ab, taxa = "Class")

tax_bar_pq(data_fungi_ab, taxa = "Class", percent_bar = TRUE)

tax_bar_pq(data_fungi_ab, taxa = "Class", fact = "Time")
#> Warning: Removed 770 rows containing missing values or values outside the scale range
#> (`geom_bar()`).

# }
```
