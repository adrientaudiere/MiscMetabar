# Test if MultiQC is installed.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_multiqc_installed(path = "multiqc")
```

## Arguments

- path:

  (default: multiqc) Path to MultiQC

## Value

A logical that say if MultiQC is installed.

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_multiqc_installed()
#> [1] TRUE
```
