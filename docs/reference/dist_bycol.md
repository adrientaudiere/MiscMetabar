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

## Examples

``` r
m1 <- matrix(runif(9), nrow = 3)
m2 <- matrix(runif(9), nrow = 3)
dist_bycol(m1, m2, nperm = 9)
#> $obs
#> [1] 0.6654949 0.2127906 0.3030498
#> 
#> $null
#> $null$length
#> [1] 0.7918465 0.4404849 0.4076450
#> 
#> $null[[2]]
#> [1] 0.7918465 0.4404849 0.4076450
#> 
#> $null[[3]]
#> [1] 0.6654949 0.2127906 0.3030498
#> 
#> $null[[4]]
#> [1] 0.7918465 0.4404849 0.4076450
#> 
#> $null[[5]]
#> [1] 0.3582721 0.4404849 0.3030498
#> 
#> $null[[6]]
#> [1] 0.3582721 0.4404849 0.3030498
#> 
#> $null[[7]]
#> [1] 0.6654949 0.2127906 0.3030498
#> 
#> $null[[8]]
#> [1] 0.3582721 0.4404849 0.3030498
#> 
#> $null[[9]]
#> [1] 0.6654949 0.2127906 0.3030498
#> 
#> 
```
