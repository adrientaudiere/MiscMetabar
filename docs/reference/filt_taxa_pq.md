# Filter taxa of a phyloseq object based on the minimum number of sequences/samples

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wraper of
[`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_pq.md).

## Usage

``` r
filt_taxa_pq(
  physeq,
  min_nb_seq = NULL,
  min_occurence = NULL,
  combination = "AND",
  clean_pq = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- min_nb_seq:

  (int default NULL) minimum number of sequences by taxa.

- min_occurence:

  (int default NULL) minimum number of sample by taxa.

- combination:

  Either "AND" (default) or "OR". If set to "AND" and both min_nb_seq
  and min_occurence are not NULL, the taxa must match the two condition
  to passe the filter. If set to "OR", taxa matching only one condition
  are kept.

- clean_pq:

  (logical) If set to TRUE, empty samples and empty taxa (ASV, OTU) are
  discarded after filtering.

## Value

a new phyloseq object

## Author

Adrien Taudière

## Examples

``` r
filt_taxa_pq(data_fungi, min_nb_seq = 20)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 32
#> Number of kept ASV 1388
#> Number of kept samples 185
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1388 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1388 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1388 reference sequences ]
filt_taxa_pq(data_fungi, min_occurence = 2)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 206
#> Number of kept ASV 1214
#> Number of kept samples 185
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1214 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1214 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1214 reference sequences ]
filt_taxa_pq(data_fungi,
  min_occurence = 2,
  min_nb_seq = 10, clean_pq = FALSE
)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 333
#> Number of kept ASV 1087
#> Number of kept samples 185
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1087 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1087 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1087 reference sequences ]
filt_taxa_pq(data_fungi,
  min_occurence = 2,
  min_nb_seq = 10,
  combination = "OR"
)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 0
#> Number of kept ASV 1420
#> Number of kept samples 185
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
```
