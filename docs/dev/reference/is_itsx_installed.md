# Test if ITSx is installed

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for testthat and examples compilation for R CMD CHECK and test
coverage.

## Usage

``` r
is_itsx_installed(path = find_itsx(), args_before_itsx = "")
```

## Arguments

- path:

  (default:
  [`find_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsx.md))
  Path to ITSx.

- args_before_itsx:

  (String, default "") A one line bash command run before ITSx, e.g. the
  value returned by
  [`install_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsx.md):
  `"source ~/miniforge3/etc/profile.d/conda.sh && conda activate itsxenv && "`.

## Value

A logical that says if ITSx can be run.

## See also

[`find_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsx.md),
[`install_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsx.md),
[`itsx_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsx_pq.md)

## Author

Adrien Taudière

## Examples

``` r
MiscMetabar::is_itsx_installed()
#> [1] FALSE
```
