# Subset samples using a conditional boolean vector.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The main objective of this function is to complete the
[`phyloseq::subset_samples()`](https://rdrr.io/pkg/phyloseq/man/subset_samples-methods.html)
function by propose a more easy (but more prone to error) way of
subset_samples. It replace the subsetting expression which used the name
of the variable in the sam_data by a boolean vector.

Warnings: you must verify the result of this function as the boolean
condition must match the order of samples in the `sam_data` slot.

This function is robust when you use the sam_data slot of the phyloseq
object used in physeq (see examples)

## Usage

``` r
subset_samples_pq(physeq, condition)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- condition:

  A boolean vector to subset samples. Length must fit the number of
  samples

## Value

a new phyloseq object

## Author

Adrien Taudière

## Examples

``` r
cond_samp <- grepl("A1", data_fungi@sam_data[["Sample_names"]])
subset_samples_pq(data_fungi, cond_samp)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 9 samples ]
#> sample_data() Sample Data:       [ 9 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]

subset_samples_pq(data_fungi, data_fungi@sam_data[["Height"]] == "Low")
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 45 samples ]
#> sample_data() Sample Data:       [ 45 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
