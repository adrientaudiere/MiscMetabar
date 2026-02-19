# Track the number of reads (= sequences), samples and cluster (e.g. ASV) for each sample

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Accept all input types supported by
[`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md):
phyloseq objects, matrices (samples x clusters), dada-class,
derep-class, lists of dada-class or derep-class, and character vectors
of fastq/fastq.gz file paths. More information are available in the
manual of the function
[`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)

## Usage

``` r
track_wkflow_samples(list_of_objects, ...)
```

## Arguments

- list_of_objects:

  (required) a list of objects passed on to
  [`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md).
  Accepts phyloseq, matrix, dada-class, derep-class, lists of dada-class
  or derep-class, and character vectors of file paths.

- ...:

  Other args passed on to
  [`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)

## Value

A list of dataframe. cf
[`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow.md)
for more information

## Author

Adrien Taudière

## Examples

``` r
tree_A10_005 <- subset_samples(data_fungi, Tree_name == "A10-005")
if (requireNamespace("pbapply")) {
  track_wkflow_samples(tree_A10_005)
}
#> You select 1 of 3 samples and conserved 67 out of 1420 taxa represented by 13314 sequences (out of 56868 sequences [23%])
#> Compute the number of sequences
#> Start object of class: phyloseq
#> Compute the number of clusters
#> Start object of class: phyloseq
#> Compute the number of samples
#> Start object of class: phyloseq
#> You select 1 of 3 samples and conserved 70 out of 1420 taxa represented by 39593 sequences (out of 56868 sequences [70%])
#> Compute the number of sequences
#> Start object of class: phyloseq
#> Compute the number of clusters
#> Start object of class: phyloseq
#> Compute the number of samples
#> Start object of class: phyloseq
#> You select 1 of 3 samples and conserved 54 out of 1420 taxa represented by 3961 sequences (out of 56868 sequences [7%])
#> Compute the number of sequences
#> Start object of class: phyloseq
#> Compute the number of clusters
#> Start object of class: phyloseq
#> Compute the number of samples
#> Start object of class: phyloseq
#> $`A10-005-B_S188_MERGED.fastq.gz`
#>   nb_sequences nb_clusters nb_samples
#> 1        13314          67          1
#> 
#> $`A10-005-H_S189_MERGED.fastq.gz`
#>   nb_sequences nb_clusters nb_samples
#> 1        39593          70          1
#> 
#> $`A10-005-M_S190_MERGED.fastq.gz`
#>   nb_sequences nb_clusters nb_samples
#> 1         3961          54          1
#> 
```
