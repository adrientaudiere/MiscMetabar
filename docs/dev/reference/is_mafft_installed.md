# Is MAFFT available?

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Check whether the MAFFT command-line aligner used by
`align_pq(method = "mafft")` is installed. Useful to guard examples,
tests and vignette chunks.

## Usage

``` r
is_mafft_installed(path = NULL)
```

## Arguments

- path:

  Optional path to the MAFFT executable. Default to NULL, in which case
  MAFFT is looked up in two places, in order: the
  `MiscMetabar.mafftpath` option, then the system `PATH`.

## Value

A logical of length one. FALSE when the `ips` package is not installed,
so that the check also guards the R-level dependency.

## See also

[`align_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/align_pq.md),
[`is_vsearch_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_vsearch_installed.md)

## Author

Adrien Taudière

## Examples

``` r
is_mafft_installed()
#> [1] TRUE
```
