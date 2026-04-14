# Install MMseqs2 from GitHub releases

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Downloads a pre-compiled MMseqs2 binary from
<https://mmseqs.com/latest/> and places it in the user data directory
for this package. Subsequent calls to
[`find_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_mmseqs2.md)
will find it automatically.

## Usage

``` r
install_mmseqs2(
  version = "latest",
  path = tools::R_user_dir("MiscMetabar", "data"),
  force = FALSE
)
```

## Arguments

- version:

  Character. Either `"latest"` (default) or a specific release tag (e.g.
  `"17-b804f"`).

- path:

  Destination directory (default:
  `tools::R_user_dir("MiscMetabar", "data")`).

- force:

  Logical. Re-download even if already installed?

## Value

The path to the installed binary (invisibly).

## See also

[`find_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_mmseqs2.md),
[`is_mmseqs2_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_mmseqs2_installed.md),
[`assign_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_mmseqs2.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
install_mmseqs2()
install_mmseqs2(force = TRUE)
} # }
```
