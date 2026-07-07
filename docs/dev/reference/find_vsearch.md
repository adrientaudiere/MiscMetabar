# Find the vsearch binary

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Searches for the vsearch binary in the following order:

1.  The `MiscMetabar.vsearchpath` option (if set)

2.  A previously installed copy in the MiscMetabar user data directory
    (via
    [`install_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_vsearch.md))

3.  The system PATH

## Usage

``` r
find_vsearch()
```

## Value

A character string with the path to the vsearch binary, or `"vsearch"`
as a fallback (relying on PATH resolution).

## See also

[`install_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_vsearch.md),
[`is_vsearch_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_vsearch_installed.md)

## Author

Adrien Taudière

## Examples

``` r
find_vsearch()
#> [1] "/home/adrien/.local/share/R/MiscMetabar/bin/vsearch"
```
