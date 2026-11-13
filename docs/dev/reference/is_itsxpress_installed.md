# Test if ITSxpress is installed

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for `testthat` and examples compilation for `R CMD check` and
test coverage.

## Usage

``` r
is_itsxpress_installed(path = find_itsxpress(), args_before_itsxpress = "")
```

## Arguments

- path:

  (default:
  [`find_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsxpress.md))
  Path to `itsxpress`.

- args_before_itsxpress:

  (String, default `""`) A one-line bash command run before `itsxpress`,
  e.g. the value returned by
  [`install_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsxpress.md):
  `"source ~/miniforge3/etc/profile.d/conda.sh && conda activate itsxpressenv && "`.

## Value

A logical indicating whether `itsxpress` can be run.

## See also

[`find_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsxpress.md),
[`install_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsxpress.md),
[`itsxpress_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsxpress_pq.md)

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_itsxpress_installed()
#> [1] FALSE
```
