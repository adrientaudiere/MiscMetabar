# Test if fastp is installed.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage.

## Usage

``` r
is_fastp_installed(path = "fastp")
```

## Arguments

- path:

  (default: fastp) Path to fastp.

## Value

A logical that say if fastp is installed.

## See also

[`fastp()`](https://adrientaudiere.github.io/MiscMetabar/reference/fastp.md)

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_fastp_installed()
#> [1] FALSE
```
