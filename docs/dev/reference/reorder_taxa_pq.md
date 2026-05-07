# Reorder taxa in otu_table/tax_table/refseq slot of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that the taxa order in a physeq object with a tree is locked by the
order of leaf in the phylogenetic tree.

## Usage

``` r
reorder_taxa_pq(physeq, names_ordered, remove_phy_tree = FALSE)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- names_ordered:

  (required) Names of the taxa (must be the same as taxa in
  `taxa_names(physeq)`) in a given order

- remove_phy_tree:

  (logical, default FALSE) If TRUE, the phylogenetic tree is removed. It
  is

## Value

A phyloseq object

## Author

Adrien Taudière

## Examples

``` r
data_fungi_ordered_by_genus <- reorder_taxa_pq(
  data_fungi,
  taxa_names(data_fungi)[order(as.vector(data_fungi@tax_table[, "Genus"]))]
)

data_fungi_mini_asc_ordered_by_abundance <- reorder_taxa_pq(
  data_fungi_mini,
  taxa_names(data_fungi_mini)[order(taxa_sums(data_fungi_mini))]
)
```
