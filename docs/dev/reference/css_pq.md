# Cumulative Sum Scaling (CSS) normalization of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wrapper around
[`metagenomeSeq::cumNorm()`](https://rdrr.io/pkg/metagenomeSeq/man/cumNorm.html)
/
[`metagenomeSeq::MRcounts()`](https://rdrr.io/pkg/metagenomeSeq/man/MRcounts.html)
implementing Cumulative Sum Scaling (Paulson et al. 2013,
[doi:10.1038/nmeth.2658](https://doi.org/10.1038/nmeth.2658) ).

## Usage

``` r
css_pq(physeq, log = TRUE)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- log:

  (logical, default `TRUE`) whether to return \\\log_2(x + 1)\\
  transformed counts (as recommended by the metagenomeSeq authors).

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with a CSS normalised `otu_table`.

## See also

[`metagenomeSeq::cumNorm()`](https://rdrr.io/pkg/metagenomeSeq/man/cumNorm.html)

## Author

Adrien Taudière

## Examples

``` r
data_f_css <- css_pq(data_fungi_mini)
#> Warning: cumNormStatFast() failed (likely samples with <=1 feature); falling back to cumNormStat(). Original error: Warning sample with one or zero features
#> Default value being used.
```
