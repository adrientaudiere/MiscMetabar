# Test if krona is installed.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_krona_installed(path = "ktImportKrona")
```

## Arguments

- path:

  (default: krona) Path to krona

## Value

A logical that say if krona is install in

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_krona_installed()
#> [1] FALSE
```
