# Test if mumu is installed.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_mumu_installed(path = "mumu")
```

## Arguments

- path:

  (default: mumu) Path to mumu

## Value

A logical that say if mumu is install in

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_mumu_installed()
#> Warning: running command 'mumu 2>&1' had status 1
#> [1] TRUE
```
