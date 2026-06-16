# Test if falco is installed.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_falco_installed(path = "falco")
```

## Arguments

- path:

  (default: falco) Path to falco

## Value

A logical that say if falco is install in

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_falco_installed()
#> [1] FALSE
```
