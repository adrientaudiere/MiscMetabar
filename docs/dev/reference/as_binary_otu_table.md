# Transform the otu_table of a [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object into a [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object with a binary otu_table.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful to test if the results are not biased by sequences bias that
appended during PCR or NGS pipeline.

## Usage

``` r
as_binary_otu_table(physeq, min_number = 1)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- min_number:

  (int) the minimum number of sequences to put a 1 in the OTU table.

## Value

A `physeq` object with only 0/1 in the OTU table

## Author

Adrien Taudière

## Examples

``` r
data(enterotype)
enterotype_bin <- as_binary_otu_table(enterotype)
```
