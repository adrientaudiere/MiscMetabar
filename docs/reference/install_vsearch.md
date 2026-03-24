# Install vsearch binary

Downloads and installs the vsearch binary from
[GitHub](https://github.com/torognes/vsearch/releases) into the
MiscMetabar user data directory. This is especially useful on Windows
where vsearch is not available from a system package manager.

After installation, all MiscMetabar functions that use vsearch will find
the binary automatically via
[`find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_vsearch.md).

## Usage

``` r
install_vsearch(
  version = "latest",
  path = tools::R_user_dir("MiscMetabar", "data"),
  force = FALSE
)
```

## Arguments

- version:

  (default: "latest") The vsearch version to install (e.g. `"2.30.5"`).
  Use `"latest"` to fetch the most recent release.

- path:

  (default: `tools::R_user_dir("MiscMetabar", "data")`) Directory where
  vsearch will be installed.

- force:

  (default: FALSE) If `TRUE`, re-download even if vsearch is already
  installed.

## Value

The path to the installed vsearch binary (invisibly).

## See also

[`find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_vsearch.md),
[`is_vsearch_installed()`](https://adrientaudiere.github.io/MiscMetabar/reference/is_vsearch_installed.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
install_vsearch()
install_vsearch(version = "2.30.5")
} # }
```
