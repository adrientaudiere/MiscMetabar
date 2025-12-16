# Filter taxa by cleaning taxa with NA at given taxonomic rank(s)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wrapper of subset_taxa_pq()

## Usage

``` r
filt_taxa_wo_NA(
  physeq,
  taxa_ranks = NULL,
  n_NA = 0,
  verbose = TRUE,
  NA_equivalent = NULL,
  clean_pq = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxa_ranks:

  A vector of taxonomic ranks. For examples c("Family","Genus"). If
  taxa_ranks is NULL (default), all ranks are used, i.e. all taxa with
  at least 1 NA will be filtered out. Numeric position of taxonomic
  ranks can also be used.

- n_NA:

  (int default = 0). Number of allowed NA by taxa in the list of the
  taxonomic ranks

- verbose:

  (logical). If TRUE, print additional information.

- NA_equivalent:

  (vector of character, default NULL). Exact matching of the character
  listed in the vector are converted as NA before to filter out taxa.

- clean_pq:

  (logical, default TRUE) If set to TRUE, empty samples are discarded
  after filtering. See
  [`clean_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/clean_pq.md).

## Value

An object of class phyloseq

## See also

[`subset_taxa_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/subset_taxa_pq.md)

## Author

Adrien Taudière

## Examples

``` r
data_fungi_wo_NA <- filt_taxa_wo_NA(data_fungi)
#> Cleaning suppress 0 taxa (  ) and 1 sample(s) ( N23-002-M_S132_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 771
#> Number of kept ASV 649
#> Number of kept samples 184
#> You filtered out 771 taxa, leading to a phyloseq object including 649 taxa without NA in the taxonomic ranks: 1 2 3 4 5 6 7 8 9 10 11 12.
filt_taxa_wo_NA(data_fungi, n_NA = 1)
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 410
#> Number of kept ASV 1010
#> Number of kept samples 185
#> You filtered out 410 taxa, leading to a phyloseq object including 1010 taxa without NA in the taxonomic ranks: 1 2 3 4 5 6 7 8 9 10 11 12.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1010 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1010 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1010 reference sequences ]
filt_taxa_wo_NA(data_fungi, taxa_ranks = c(1:3))
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 93
#> Number of kept ASV 1327
#> Number of kept samples 185
#> You filtered out 93 taxa, leading to a phyloseq object including 1327 taxa without NA in the taxonomic ranks: 1 2 3.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1327 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1327 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1327 reference sequences ]

filt_taxa_wo_NA(data_fungi, taxa_ranks = c("Trait", "Confidence.Ranking"))
#> Cleaning suppress 0 taxa (  ) and 0 sample(s) (  ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 0
#> Number of kept ASV 1420
#> Number of kept samples 185
#> You filtered out 0 taxa, leading to a phyloseq object including 1420 taxa without NA in the taxonomic ranks: Trait Confidence.Ranking.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
filt_taxa_wo_NA(data_fungi,
  taxa_ranks = c("Trait", "Confidence.Ranking"),
  NA_equivalent = c("-", "NULL")
)
#> Cleaning suppress 0 taxa (  ) and 4 sample(s) ( N23-002-M_S132_MERGED.fastq.gz / O26-004-M_S150_MERGED.fastq.gz / Y28-002-H_S179_MERGED.fastq.gz / Y28-002-M_S180_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 1420
#> Number of filtered-out ASV 1198
#> Number of kept ASV 222
#> Number of kept samples 181
#> You filtered out 1198 taxa, leading to a phyloseq object including 222 taxa without NA in the taxonomic ranks: Trait Confidence.Ranking.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 222 taxa and 181 samples ]
#> sample_data() Sample Data:       [ 181 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 222 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 222 reference sequences ]
```
