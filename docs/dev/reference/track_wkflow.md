# Track the number of reads (= sequences), samples and cluster (e.g. ASV) from various objects including dada-class and derep-class.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

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
  compute_occurrences = FALSE,
  compute_taxo_info = FALSE,
  ranks_for_na_prop = c("Class", "Order", "Family", "Genus"),
  compute_sam_metadata = FALSE,
  compute_seq_length = FALSE,
  compute_genetic_diversity = FALSE,
  compute_factor_counts = FALSE,
  factor = NULL,
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

- compute_occurrences:

  (logical, default FALSE) If TRUE, add a column `nb_occurrences` with
  the total number of non-zero cells in the OTU table (number of (taxon,
  sample) presences), computed via
  [`as_binary_otu_table()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/as_binary_otu_table.md).
  Phyloseq objects only (NA otherwise).

- compute_taxo_info:

  (logical, default FALSE) If TRUE, add a column `nb_rank` (number of
  taxonomic ranks,
  [`phyloseq::rank_names()`](https://rdrr.io/pkg/phyloseq/man/rank_names.html))
  and one column `prop_na_<rank>` per rank in `ranks_for_na_prop` giving
  the proportion of NA assignments among present taxa.

- ranks_for_na_prop:

  (character, default `c("Class","Order","Family","Genus")`) Rank names
  for which the proportion of NA is computed when
  `compute_taxo_info = TRUE`. Ignored otherwise.

- compute_sam_metadata:

  (logical, default FALSE) If TRUE, add a column `nb_sam_metadata` with
  the number of columns in `sample_data`.

- compute_seq_length:

  (logical, default FALSE) If TRUE, add three columns `mean_length_seq`,
  `max_length_seq`, `min_length_seq` from the width of the reference
  sequences (`Biostrings::width()`) of present taxa.

- compute_genetic_diversity:

  (logical, default FALSE) If TRUE, add two columns
  `genetic_diversity_weighted` (abundance-weighted nucleotide diversity,
  pi) and `genetic_diversity_unweighted` (mean pairwise distance), both
  computed from `@ref_seq` of present taxa with
  [`DECIPHER::DistanceMatrix()`](https://rdrr.io/pkg/DECIPHER/man/DistanceMatrix.html),
  which handles unequal-length sequences by using the shorter length in
  each pairwise comparison. Requires the DECIPHER package (in
  `Suggests`). May be slow for objects with many taxa.

- compute_factor_counts:

  (logical, default FALSE) If TRUE, add one column `n_samples_<level>`
  per level of `factor` with the number of samples in each level.

- factor:

  (character, default NULL) Name of a variable in `sample_data` used
  when `compute_factor_counts = TRUE`.

- ...:

  Additional arguments passed on to
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/clean_pq.md)
  function.

## Value

The number of sequences, clusters (e.g. OTUs, ASVs) and samples for each
object. When any `compute_*` extra-metrics argument is TRUE, extra
columns are appended (phyloseq objects only; NA for other object types).
These are propagated to
[`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/track_wkflow_samples.md)
via `...`.

## See also

[`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/track_wkflow_samples.md)

## Author

Adrien Taudière

## Examples

``` r
data(enterotype)
if (requireNamespace("pbapply")) {
  track_wkflow(list(data_fungi_mini, enterotype), taxonomy_rank = c(3, 5))
  track_wkflow(list(
    "data FUNGI" = data_fungi_mini,
    "fastq files forward" =
      unlist(list_fastq_files(system.file("extdata", package = "MiscMetabar"),
        paired_end = FALSE
      ))
  ))
  track_wkflow(list(data_fungi_mini, enterotype),
    compute_occurrences = TRUE, compute_taxo_info = TRUE,
    compute_seq_length = TRUE, compute_genetic_diversity = TRUE)
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
#> Compute the number of sequences
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Compute the number of clusters
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#> Compute the number of samples
#> Start object of class: phyloseq
#> Start object of class: phyloseq
#>   nb_sequences nb_clusters nb_samples nb_occurrences nb_rank prop_na_Class
#> 1  569525.0000          45        137            560      12    0.02222222
#> 2     279.9807         553        280              0       1            NA
#>   prop_na_Order prop_na_Family prop_na_Genus mean_length_seq max_length_seq
#> 1    0.06666667      0.1333333   0.177777778        343.0889            387
#> 2            NA             NA   0.003616637              NA             NA
#>   min_length_seq genetic_diversity_weighted genetic_diversity_unweighted
#> 1            303                  0.4743904                    0.4924878
#> 2             NA                         NA                           NA
```
