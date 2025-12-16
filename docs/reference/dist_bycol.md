# Compute paired distances among matrix (e.g. otu_table)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

May be used to verify ecological distance among samples.

## Usage

``` r
dist_bycol(x, y, method = "bray", nperm = 99, ...)
```

## Arguments

- x:

  (required) A first matrix.

- y:

  (required) A second matrix.

- method:

  (default: 'bray') the method to use internally in the vegdist
  function.

- nperm:

  (int) The number of permutations to perform.

- ...:

  Additional arguments passed on
  to[`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function

## Value

A list of length two : (i) a vector of observed distance (\$obs) and
(ii) a matrix of the distance after randomization (\$null)

## Note

the first column of the first matrix is compare to the first column of
the second matrix, the second column of the first matrix is compare to
the second column of the second matrix and so on.

## See also

[`vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)

## Author

Adrien Taudière
