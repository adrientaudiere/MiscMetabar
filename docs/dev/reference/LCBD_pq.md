# Compute and test local contributions to beta diversity (LCBD) of samples

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
function in the case of `physeq` object.

## Usage

``` r
LCBD_pq(physeq, p_adjust_method = "BH", ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- p_adjust_method:

  (chr, default "BH"): the method used to adjust p-value

- ...:

  Additional arguments passed on to
  [`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
  function

## Value

An object of class `beta.div` see
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
function for more information

## See also

[plot_LCBD_pq](https://adrientaudiere.github.io/MiscMetabar/dev/reference/plot_LCBD_pq.md),
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)

## Author

Adrien Taudière This function is mainly a wrapper of the work of others.
Please make a reference to
[`adespatial::beta.div()`](http://adeverse.github.io/adespatial/reference/beta.div.md)
if you use this function.

## Examples

``` r
# \donttest{
if (requireNamespace("adespatial")) {
  data_f <- clean_pq(prune_samples(
    sample_names(data_fungi_sp_known)[1:10],
    data_fungi_sp_known
  ))
  res <- LCBD_pq(data_f, nperm = 5)
  str(res)
  length(res$LCBD)
  length(res$SCBD)
}
#> Cleaning suppress 412 taxa and 0 samples.
#> List of 8
#>  $ beta  : Named num [1:2] 8.267 0.919
#>   ..- attr(*, "names")= chr [1:2] "SStotal" "BDtotal"
#>  $ SCBD  : Named num [1:239] 0.0123 0.062 0.0302 0.0438 0.0288 ...
#>   ..- attr(*, "names")= chr [1:239] "ASV2" "ASV8" "ASV12" "ASV18" ...
#>  $ LCBD  : Named num [1:10] 0.0833 0.1045 0.0866 0.1044 0.0839 ...
#>   ..- attr(*, "names")= chr [1:10] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  $ p.LCBD: Named num [1:10] 0.833 0.667 0.833 0.833 1 ...
#>   ..- attr(*, "names")= chr [1:10] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  $ p.adj : Named num [1:10] 0.926 0.926 0.926 0.926 1 ...
#>   ..- attr(*, "names")= chr [1:10] "A10-005-B_S188_MERGED.fastq.gz" "A10-005-H_S189_MERGED.fastq.gz" "A10-005-M_S190_MERGED.fastq.gz" "A12-007_S191_MERGED.fastq.gz" ...
#>  $ method: chr [1:2] "hellinger" NA
#>  $ note  : chr "Info -- This coefficient is Euclidean"
#>  $ D     : logi NA
#>  - attr(*, "class")= chr "beta.div"
#> [1] 239
# }
```
