# Check whether MMseqs2 is installed and callable

Tries to run `mmseqs version` and returns `TRUE` if it succeeds.

## Usage

``` r
is_mmseqs2_installed(path = find_mmseqs2())
```

## Arguments

- path:

  Path to the `mmseqs` binary (default:
  [`find_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_mmseqs2.md)).

## Value

Logical.

## See also

[`find_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_mmseqs2.md),
[`install_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/install_mmseqs2.md),
[`assign_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_mmseqs2.md)

## Author

Adrien Taudière

## Examples

``` r
is_mmseqs2_installed()
#> [1] TRUE
```
