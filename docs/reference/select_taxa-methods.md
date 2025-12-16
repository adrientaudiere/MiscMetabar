# Select a subset of taxa in a specified order where possible

Select (a subset of) taxa; if `x` allows taxa to be reordered, then taxa
are given in the specified order.

## Usage

``` r
select_taxa(x, taxa, reorder = TRUE)

# S4 method for class 'sample_data,character'
select_taxa(x, taxa)

# S4 method for class 'otu_table,character'
select_taxa(x, taxa, reorder = TRUE)

# S4 method for class 'taxonomyTable,character'
select_taxa(x, taxa, reorder = TRUE)

# S4 method for class 'XStringSet,character'
select_taxa(x, taxa, reorder = TRUE)

# S4 method for class 'phylo,character'
select_taxa(x, taxa)

# S4 method for class 'phyloseq,character'
select_taxa(x, taxa, reorder = TRUE)
```

## Arguments

- x:

  A phyloseq object or phyloseq component object

- taxa:

  Character vector of taxa to select, in requested order

- reorder:

  Logical specifying whether to use the order in `taxa` (TRUE) or keep
  the order in `taxa_names(x)` (FALSE)

## Details

This is a simple selector function that is like `prune_taxa(taxa, x)`
when `taxa` is a character vector but always gives the taxa in the order
`taxa` if possible (that is, except for phy_tree's and phyloseq objects
that contain phy_tree's).

## Author

Michael R. McLaren (orcid:
[0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
