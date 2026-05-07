# Add dna in `refseq` slot of a `physeq` object using taxa names and renames taxa using prefix_taxa_names and number (default Taxa_1, Taxa_2 ...)

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful in targets bioinformatic pipeline.

## Usage

``` r
add_dna_to_phyloseq(physeq, prefix_taxa_names = "Taxa_")
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- prefix_taxa_names:

  (default "Taxa\_"): the prefix of taxa names (eg. "ASV\_" or "OTU\_")

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with `refseq` slot and new taxa names

## Author

Adrien Taudière

## Examples

``` r
pq_seq_names <- phyloseq::phyloseq(
  phyloseq::otu_table(data_fungi_mini),
  phyloseq::sample_data(data_fungi_mini),
  phyloseq::tax_table(data_fungi_mini)
)
phyloseq::taxa_names(pq_seq_names) <- as.character(phyloseq::refseq(data_fungi_mini))
add_dna_to_phyloseq(pq_seq_names)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]
```
