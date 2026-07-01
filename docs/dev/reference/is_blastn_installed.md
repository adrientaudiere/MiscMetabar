# Test if blastn is installed.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_blastn_installed(path = "blastn")
```

## Arguments

- path:

  (default: blastn) Path to blastn (NCBI BLAST+)

## Value

A logical that say if blastn is installed.

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_blastn_installed()
#> [1] TRUE
```
