# Test if swarm is installed.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage

## Usage

``` r
is_swarm_installed(path = "swarm")
```

## Arguments

- path:

  (default: swarm) Path to falco

## Value

A logical that say if swarm is install in

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_swarm_installed()
#> [1] TRUE
```
