# Find the ITSx executable

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Resolution order: the `MiscMetabar.itsxpath` option (if set), then the
system `PATH`. ITSx is usually installed in a conda environment (see
[`install_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsx.md));
in that case it is not on the `PATH` of R and the environment is
activated with the `args_before_itsx` argument of
[`itsx_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsx_pq.md)
and
[`is_itsx_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsx_installed.md)
instead.

## Usage

``` r
find_itsx()
```

## Value

A character string with the path to ITSx, or `"ITSx"` as a fallback
(relying on `PATH` resolution, possibly after `args_before_itsx`).

## See also

[`itsx_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsx_pq.md),
[`is_itsx_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsx_installed.md),
[`install_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsx.md)

## Author

Adrien Taudière

## Examples

``` r
find_itsx()
#> [1] "ITSx"
```
