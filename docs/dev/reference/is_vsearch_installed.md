# Test if vsearch is installed.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_vsearch_installed(path = find_vsearch())
```

## Arguments

- path:

  (default:
  [`find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_vsearch.md))
  Path to vsearch

## Value

A logical that say if vsearch is install in

## See also

[`find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_vsearch.md),
[`install_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_vsearch.md)

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_vsearch_installed()
#> [1] TRUE
```
