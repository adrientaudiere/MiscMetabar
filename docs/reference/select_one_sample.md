# Select one sample from a physeq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mostly for internal used, for example in function
[`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow_samples.md).

## Usage

``` r
select_one_sample(physeq, sam_name, silent = FALSE)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- sam_name:

  (required) The sample name to select

- silent:

  (logical) If true, no message are printing.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with one sample

## Author

Adrien Taudière

## Examples

``` r
A8_005 <- select_one_sample(data_fungi, "A8-005_S4_MERGED.fastq.gz")
#> You select 1 of 185 samples and conserved 83 out of 1420 taxa represented by 13875 sequences (out of 1839124 sequences [1%])
A8_005
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 83 taxa and 1 samples ]
#> sample_data() Sample Data:       [ 1 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 83 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 83 reference sequences ]
```
