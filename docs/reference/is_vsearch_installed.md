# Test if vsearch is installed.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_vsearch_installed(path = "vsearch")
```

## Arguments

- path:

  (default: vsearch) Path to vsearch

## Value

A logical that say if vsearch is install in

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_vsearch_installed()
#> [1] TRUE
```
