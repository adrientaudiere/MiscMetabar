# Load a phyloseq object from an RData file written by `save_pq()` or `write_pq()`

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This is the reverse function of
[`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/save_pq.md)
and
[`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/write_pq.md)
when used with `rdata = TRUE`. Those functions save the phyloseq object
in an RData file named `physeq.RData` (the object is always stored under
the name `physeq`). Contrary to
[`base::load()`](https://rdrr.io/r/base/load.html), `load_pq()` returns
the phyloseq object so it can be assigned to any name, without modifying
the calling environment.

## Usage

``` r
load_pq(path = NULL)
```

## Arguments

- path:

  (required) a path to the folder where the phyloseq object was saved,
  or directly the path to the `.RData`/`.rda` file. If a folder is
  given, the file `physeq.RData` inside it is loaded.

## Value

The loaded phyloseq object.

## See also

[`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/save_pq.md),
[`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/write_pq.md),
[`read_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/read_pq.md)

## Author

Adrien Taudière

## Examples

``` r
save_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"))
# Load the object and assign it to any name you like
my_pq <- load_pq(path = paste0(tempdir(), "/phyloseq"))
unlink(paste0(tempdir(), "/phyloseq"), recursive = TRUE)
```
