# iNterpolation and EXTrapolation of Hill numbers (with iNEXT)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that this function is quite time-consuming due to high
dimensionality in metabarcoding community matrix.

## Usage

``` r
iNEXT_pq(physeq, merge_sample_by = NULL, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- merge_sample_by:

  (default: NULL) if not `NULL` samples of physeq are merged using the
  vector set by `merge_sample_by`. This merging used the
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md).
  In the case of
  [`biplot_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/biplot_pq.md)
  this must be a factor with two levels only.

- ...:

  Other arguments for the
  [`iNEXT::iNEXT()`](https://rdrr.io/pkg/iNEXT/man/iNEXT.html) function

## Value

see [`iNEXT::iNEXT()`](https://rdrr.io/pkg/iNEXT/man/iNEXT.html)
documentation

## Author

Adrien Taudière This function is mainly a wrapper of the work of others.
Please make a reference to
[`iNEXT::iNEXT()`](https://rdrr.io/pkg/iNEXT/man/iNEXT.html) if you use
this function.

## Examples

``` r
# \donttest{
if (requireNamespace("iNEXT")) {
  data("GlobalPatterns", package = "phyloseq")
  GPsubset <- subset_taxa(
    GlobalPatterns,
    GlobalPatterns@tax_table[, 1] == "Bacteria"
  )
  GPsubset <- subset_taxa(
    GPsubset,
    rowSums(GPsubset@otu_table) > 20000
  )
  GPsubset <- subset_taxa(
    GPsubset,
    rowSums(is.na(GPsubset@tax_table)) == 0
  )
  GPsubset@sam_data$human <- GPsubset@sam_data$SampleType %in%
    c("Skin", "Feces", "Tong")
  res_iNEXT <- iNEXT_pq(
    GPsubset,
    merge_sample_by = "human",
    q = 1,
    datatype = "abundance",
    nboot = 2
  )
  iNEXT::ggiNEXT(res_iNEXT)
  # iNEXT::ggiNEXT(res_iNEXT, type = 2)
  # iNEXT::ggiNEXT(res_iNEXT, type = 3)
}
#> Error in eval(e, x, parent.frame()): object 'GPsubset' not found
# }
```
