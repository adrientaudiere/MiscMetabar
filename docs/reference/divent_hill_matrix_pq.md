# Compute Hill diversity numbers for all samples in an OTU table

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Iterates over all samples in an OTU table and computes Hill diversity
numbers using
[`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html).

## Usage

``` r
divent_hill_matrix_pq(comm, q, ...)
```

## Arguments

- comm:

  (data.frame or matrix) OTU table with samples as rows and taxa as
  columns.

- q:

  (numeric vector) Hill diversity orders to compute. Hill numbers are
  more appropriate in DNA metabarcoding studies when `q > 0` (Alberdi &
  Gilbert, 2019; Calderón-Sanou et al., 2019).

- ...:

  Additional arguments passed to
  [`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html)
  (e.g. `estimator = "naive"` to reproduce vegan-equivalent results).

## Value

A data.frame with one row per sample and one column per value in `q`.
Column names are the string representation of the `q` values. Row names
match the input row names.

## References

Alberdi, A., & Gilbert, M. T. P. (2019). A guide to the application of
Hill numbers to DNA-based diversity analyses. *Molecular Ecology
Resources*.
[doi:10.1111/1755-0998.13014](https://doi.org/10.1111/1755-0998.13014)

Calderón-Sanou, I., Münkemüller, T., Boyer, F., Zinger, L., & Thuiller,
W. (2019). From environmental DNA sequences to ecological conclusions:
How strong is the influence of methodological choices? *Journal of
Biogeography*, 47.
[doi:10.1111/jbi.13681](https://doi.org/10.1111/jbi.13681)

## See also

[`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html),
[`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md),
[`hill_tuckey_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_tuckey_pq.md)

## Examples

``` r
data("data_fungi_mini", package = "MiscMetabar")
data_f <- prune_samples(
  sample_names(data_fungi_mini)[1:5],
  data_fungi_mini
)
otu <- as.data.frame(phyloseq::otu_table(
  taxa_as_columns(data_f)
))
#> Taxa are now in columns.
divent_hill_matrix_pq(otu, q = c(0, 1, 2))
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#> ! Function arguments cannot be checked because the package divent is not attached.
#> → Add `CheckArguments=FALSE` to suppress this warning or run `library('divent')`.
#>                                 0        1        2
#> A10-005-B_S188_MERGED.fastq.gz  6 3.778681 3.321022
#> A10-005-H_S189_MERGED.fastq.gz  7 1.012331 1.002750
#> A10-005-M_S190_MERGED.fastq.gz  7 4.151787 3.664966
#> A12-007_S191_MERGED.fastq.gz   10 1.658119 1.226407
#> A12-007-B_S2_MERGED.fastq.gz    2 1.987268 1.974830
```
