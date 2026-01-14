# Merge groups of elements within a vector by a function

Internal function used in
[`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
to merge variables. Note, owing to the use of
[`split()`](https://rdrr.io/r/base/split.html), the merged elements in
the new vector will be reordered according to `group`.

## Usage

``` r
merge_groups(x, group, f = unique_or_na)
```

## Arguments

- x:

  A vector whose elements will be merged.

- group:

  A vector such that `as.factor(group)` defines the grouping.

- f:

  A function that, when applied to a subvector of x, returns a single
  value. Can also be a formula as interpretted by
  [`purrr::as_mapper()`](https://purrr.tidyverse.org/reference/as_mapper.html).

## Author

Michael R. McLaren (orcid:
[0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
