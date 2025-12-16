# Track the number of reads (= sequences), samples and cluster (e.g. ASV) from various objects including dada-class and derep-class.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

- List of fastq and fastg.gz files -\> nb of reads and samples

- List of dada-class -\> nb of reads, clusters (ASV) and samples

- List of derep-class -\> nb of reads, clusters (unique sequences) and
  samples

- Matrix of samples x clusters (e.g. `otu_table`) -\> nb of reads,
  clusters and samples

- Phyloseq-class -\> nb of reads, clusters and samples

## Usage

``` r
track_wkflow(
  list_of_objects,
  obj_names = NULL,
  clean_pq = FALSE,
  taxonomy_rank = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- list_of_objects:

  (required) a list of objects

- obj_names:

  A list of names corresponding to the list of objects

- clean_pq:

  (logical) If set to TRUE, empty samples and empty ASV are discarded
  before clustering.

- taxonomy_rank:

  A vector of int. Define the column number of taxonomic rank
  `in physeq@tax_table` to compute the number of unique value. Default
  is NULL and do not compute values for any taxonomic rank

- verbose:

  (logical) If true, print some additional messages.

- ...:

  Additional arguments passed on to
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md)
  function.

## Value

The number of sequences, clusters (e.g. OTUs, ASVs) and samples for each
object.

## See also

[`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/reference/track_wkflow_samples.md)

## Author

Adrien Taudière

## Examples

``` r
data(enterotype)
if (requireNamespace("pbapply")) {
  track_wkflow(list(data_fungi, enterotype), taxonomy_rank = c(3, 5))
  track_wkflow(list(
    "data FUNGI" = data_fungi,
    "fastq files forward" =
      unlist(list_fastq_files(system.file("extdata", package = "MiscMetabar"),
        paired_end = FALSE
      ))
  ))
}
#> Compute the number of sequences
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Compute the number of clusters
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Compute the number of samples
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Compute the number of values in taxonomic rank
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Compute the number of sequences
#> Start object of class: phyloseq
#> Start object of class: character
#> Compute the number of clusters
#> Start object of class: phyloseq
#> Start object of class: character
#> Compute the number of samples
#> Start object of class: phyloseq
#> Start object of class: character
#>                     nb_sequences nb_clusters nb_samples
#> data FUNGI               1839124        1420        185
#> fastq files forward          100          NA          1
```
