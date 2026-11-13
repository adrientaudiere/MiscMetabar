# Install ITSxpress in a conda environment

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

ITSxpress (Rivers et al. 2018; Einarsson & Rivers 2024) is a Python tool
that requires VSEARCH, HMMER 3 and Biopython. The simplest way to
install all of them is the bioconda package, which this function
installs in a dedicated conda environment. It returns the string to pass
as `args_before_itsxpress` to
[`itsxpress_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsxpress_pq.md)
and
[`is_itsxpress_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsxpress_installed.md),
which activates that environment before each call.

## Usage

``` r
install_itsxpress(env_name = "itsxpressenv", conda = NULL, force = FALSE)
```

## Arguments

- env_name:

  (default: `"itsxpressenv"`) Name of the conda environment.

- conda:

  (default: `NULL`) Path to `mamba` or `conda`. If `NULL`, `mamba` then
  `conda` are looked up on the `PATH`.

- force:

  (default: `FALSE`) If `TRUE`, the environment is created even if one
  of that name already exists (conda then replaces it).

## Value

The `args_before_itsxpress` string activating the environment
(invisibly).

## See also

[`itsxpress_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsxpress_pq.md),
[`is_itsxpress_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsxpress_installed.md),
[`find_itsxpress()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsxpress.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
prelude <- install_itsxpress()
is_itsxpress_installed(args_before_itsxpress = prelude)
} # }
```
