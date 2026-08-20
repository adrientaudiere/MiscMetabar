# Plot the quality of reads in a fastq file

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wrapper of
[`Rsearch::plot_read_quality()`](https://cassandrahjo.github.io/Rsearch/reference/plot_read_quality.html).
Plot the per-base quality score distribution of the reads of a fastq
file, either as raw Phred scores or converted to expected-error (EE)
rates.

## Usage

``` r
plot_read_quality(
  fastq_input,
  use_ee_rate = FALSE,
  plot_title = TRUE,
  alpha = 0.5
)
```

## Arguments

- fastq_input:

  (required) Path to a fastq file (optionally gzipped).

- use_ee_rate:

  (logical, default FALSE) Convert quality scores to expected error
  rates.

- plot_title:

  (logical, default TRUE) Add a title to the plot.

- alpha:

  (numeric, default 0.5) Transparency of the quality lines.

## Value

A ggplot2-based plot object.

## See also

[`Rsearch::plot_read_quality()`](https://cassandrahjo.github.io/Rsearch/reference/plot_read_quality.html),
[`plot_ee_rate_dist()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/plot_ee_rate_dist.md)

## Author

Adrien Taudière

## Examples

``` r
plot_read_quality(
  system.file("extdata", "ex.fastq", package = "MiscMetabar")
)
```
