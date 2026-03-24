# Find the MMseqs2 binary

Looks for the MMseqs2 binary in three places, in order:

1.  The option `MiscMetabar.mmseqs2path` (if set).

2.  A local copy installed by
    [`install_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/install_mmseqs2.md)
    in the user data directory.

3.  The system `PATH`.

## Usage

``` r
find_mmseqs2()
```

## Value

A character string with the path to the `mmseqs` binary.

## See also

[`install_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/install_mmseqs2.md),
[`is_mmseqs2_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_mmseqs2_installed.md),
[`assign_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_mmseqs2.md)

## Author

Adrien Taudière
