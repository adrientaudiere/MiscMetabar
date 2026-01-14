# Get the unique value in x or NA if none

If `unique(x)` is a single value, return it; otherwise, return an NA of
the same type as `x`. If `x` is a factor, then the levels and ordered
status will be kept in either case. If `x` is a non-atomic vector (i.e.
a list), then the logical `NA` will be used.

## Usage

``` r
unique_or_na(x)
```

## Arguments

- x:

  A vector

## Value

Either a single value (if `unique(x)` return a single value) or a NA

## Author

Michael R. McLaren (orcid:
[0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))

## Examples

``` r
f <- factor(c("a", "a", "b", "c"), ordered = TRUE)
unique_or_na(f)
#> [1] <NA>
#> Levels: a < b < c
unique_or_na(f[1:2])
#> [1] a
#> Levels: a < b < c

x <- c("a", "b", "a")
unique_or_na(x[c(1, 3)])
#> [1] "a"
unique_or_na(x)
#> [1] NA
unique_or_na(x) %>% typeof()
#> [1] "character"
```
