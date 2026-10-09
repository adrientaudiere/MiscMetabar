# Install ITSx in a conda environment

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

ITSx (Bengtsson-Palme et al. 2013) is a Perl program that needs HMMER 3.
The simplest way to get both is the bioconda package, which this
function installs in a dedicated conda environment. It returns the
string to pass as `args_before_itsx` to
[`itsx_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsx_pq.md)
and
[`is_itsx_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsx_installed.md),
which activates that environment before each call.

## Usage

``` r
install_itsx(env_name = "itsxenv", conda = NULL, force = FALSE)
```

## Arguments

- env_name:

  (default: "itsxenv") Name of the conda environment.

- conda:

  (default: NULL) Path to `mamba` or `conda`. If NULL, `mamba` then
  `conda` are looked up on the `PATH`.

- force:

  (default: FALSE) If TRUE, the environment is created even if one of
  that name already exists (conda then replaces it).

## Value

The `args_before_itsx` string activating the environment (invisibly).

## Without conda

Download ITSx from <https://microbiology.se/software/itsx/>, install
HMMER 3 (<http://hmmer.org/>; `sudo apt install hmmer` on Debian /
Ubuntu, `brew install hmmer` on macOS), put the `ITSx` script and its
`ITSx_db` folder on the `PATH`, or set
`options(MiscMetabar.itsxpath = "/path/to/ITSx")`.

## See also

[`itsx_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/itsx_pq.md),
[`is_itsx_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_itsx_installed.md),
[`find_itsx()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_itsx.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
prelude <- install_itsx()
is_itsx_installed(args_before_itsx = prelude)
} # }
```
