# Plot the distribution of expected error (EE) rates of reads

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wrapper of
[`Rsearch::plot_ee_rate_dist()`](https://cassandrahjo.github.io/Rsearch/reference/plot_ee_rate_dist.html).
Plot the distribution of the expected error (EE) rate of the reads of a
fastq file, a key diagnostic before quality filtering (e.g. with `maxEE`
thresholds).

## Usage

``` r
plot_ee_rate_dist(
  fastq_input,
  n_bins = 30,
  plot_title = "Distribution of the expected error (EE) rate of reads"
)
```

## Arguments

- fastq_input:

  (required) Path to a fastq file (optionally gzipped).

- n_bins:

  (integer, default 30) Number of histogram bins.

- plot_title:

  (character) Title of the plot.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## See also

[`Rsearch::plot_ee_rate_dist()`](https://cassandrahjo.github.io/Rsearch/reference/plot_ee_rate_dist.html),
[`plot_read_quality()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/plot_read_quality.md)

## Author

Adrien Taudière

## Examples

``` r
plot_ee_rate_dist(
  system.file("extdata", "ex.fastq", package = "MiscMetabar")
)
```
