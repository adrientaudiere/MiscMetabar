# Summarize a tax_table (taxonomic slot of phyloseq object) using gtsummary

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mainly a wrapper for the
[`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)
function in the case of `physeq` object.

## Usage

``` r
tbl_sum_taxtable(physeq, taxonomic_ranks = NULL, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxonomic_ranks:

  A list of taxonomic ranks we want to summarized.

- ...:

  Additional arguments passed on to
  [`gtsummary::tbl_summary()`](https://www.danieldsjoberg.com/gtsummary/reference/tbl_summary.html)

## Value

A table of class c('tbl_summary', 'gtsummary')

## Author

Adrien Taudière

## Examples

``` r
tbl_sum_taxtable(data_fungi_mini)


  

Characteristic
```
