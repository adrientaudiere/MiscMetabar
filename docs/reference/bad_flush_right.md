# Replace all values with NA upon seeing a bad value

Helper for
[`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)

## Usage

``` r
bad_flush_right(x, bad = "BAD", na_bad = FALSE, k = length(x))
```

## Arguments

- x:

  a vector

- bad:

  the string representing a bad value

- na_bad:

  whether NAs should also be treated as bad

- k:

  the index to which values are flushed

## Author

Michael R. McLaren (orcid:
[0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
