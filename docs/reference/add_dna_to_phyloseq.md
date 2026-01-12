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
