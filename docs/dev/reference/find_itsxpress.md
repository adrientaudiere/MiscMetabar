# Find the ITSxpress executable

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Resolution order: the `MiscMetabar.itsxpresspath` option (if set), then
the system `PATH`. ITSxpress is usually installed in a conda environment
(see
[`install_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsxpress.md));
in that case it is not on the `PATH` of R and the environment is
activated with the `args_before_itsxpress` argument of
[`itsxpress_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsxpress_pq.md)
and
[`is_itsxpress_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsxpress_installed.md)
instead.

## Usage

``` r
find_itsxpress()
```

## Value

A character string with the path to `itsxpress`, or `"itsxpress"` as a
fallback (relying on `PATH` resolution, possibly after
`args_before_itsxpress`).

## See also

[`itsxpress_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsxpress_pq.md),
[`is_itsxpress_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsxpress_installed.md),
[`install_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_itsxpress.md)

## Author

Adrien Taudière

## Examples

``` r
find_itsxpress()
#> [1] "itsxpress"
```
