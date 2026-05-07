# Plot kmer complexity of references sequences of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wrapper of
[`dada2::seqComplexity()`](https://rdrr.io/pkg/dada2/man/seqComplexity.html)

## Usage

``` r
plot_complexity_pq(
  physeq,
  kmer_size = 2,
  window = NULL,
  by = 5,
  bins = 100,
  aggregate = FALSE,
  vline_random_kmer = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- kmer_size:

  int (default 2) The size of the kmers (or "oligonucleotides" or
  "words") to use.

- window:

  (int, default NULL) The width in nucleotides of the moving window. If
  NULL the whole sequence is used.

- by:

  (int, default 5) The step size in nucleotides between each moving
  window tested.

- bins:

  (int, default 100). The number of bins to use for the histogram.

- aggregate:

  (logical, default FALSE) If TRUE, compute an aggregate quality profile
  for all samples

- vline_random_kmer:

  (logical, default TRUE) If TRUE, add a vertical line at the value for
  random kmer (equal to 4^kmerSize))

- ...:

  Arguments passed on to geom_histogram.

## Value

A ggplot2 object

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`dada2::seqComplexity()`](https://rdrr.io/pkg/dada2/man/seqComplexity.html)

## See also

[`dada2::seqComplexity()`](https://rdrr.io/pkg/dada2/man/seqComplexity.html),
[`dada2::plotComplexity()`](https://rdrr.io/pkg/dada2/man/plotComplexity.html)

## Author

Adrien Taudière

## Examples

``` r
plot_complexity_pq(subset_samples(data_fungi_mini, Height == "High"),
  vline_random_kmer = FALSE
)

# plot_complexity_pq(subset_samples(data_fungi_mini, Height == "Low"),
#  aggregate = FALSE, kmer_size = 4
# )
# plot_complexity_pq(subset_samples(data_fungi, Height == "Low"),
#  kmer_size = 4)
```
