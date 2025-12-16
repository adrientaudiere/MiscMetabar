# Create sample data without adjusting row/sample names

[`phyloseq::sample_data()`](https://rdrr.io/pkg/phyloseq/man/sample_data-methods.html)
will change the sample names from the row names if they are
`as.character(seq(1, row(object)))`. This function instead keeps the
names as is.

## Usage

``` r
sample_data_stable(object)
```

## Arguments

- object:

  A "data.frame"-class object

## Author

Michael R. McLaren (orcid:
[0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))

## Examples

``` r
x <- data.frame(var1 = letters[1:3], var2 = 7:9)
rownames(x)
#> [1] "1" "2" "3"
sample_data(x)
#>     var1 var2
#> sa1    a    7
#> sa2    b    8
#> sa3    c    9
MiscMetabar:::sample_data_stable(x)
#>   var1 var2
#> 1    a    7
#> 2    b    8
#> 3    c    9
```
