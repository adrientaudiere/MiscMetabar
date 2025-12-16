# Heat tree from `metacoder` package using `tax_table` slot

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that the number of ASV is store under the name `n_obs` and the
number of sequences under the name `nb_sequences`

## Usage

``` r
heat_tree_pq(physeq, taxonomic_level = NULL, ...)
```

## Arguments

- physeq:

  (required): a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxonomic_level:

  (default: NULL): a vector of selected taxonomic level using their
  column numbers (e.g. taxonomic_level = 1:7)

- ...:

  Arguments passed on to
  [`heat_tree`](https://rdrr.io/pkg/metacoder/man/heat_tree.html)

## Value

A plot

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("metacoder")) {
  library("metacoder")
  data("GlobalPatterns", package = "phyloseq")

  GPsubset <- subset_taxa(
    GlobalPatterns,
    GlobalPatterns@tax_table[, 1] == "Bacteria"
  )

  GPsubset <- subset_taxa(
    GPsubset,
    rowSums(GPsubset@otu_table) > 5000
  )

  GPsubset <- subset_taxa(
    GPsubset,
    rowSums(is.na(GPsubset@tax_table)) == 0
  )

  heat_tree_pq(GPsubset,
    node_size = n_obs,
    node_color = n_obs,
    node_label = taxon_names,
    tree_label = taxon_names,
    node_size_trans = "log10 area"
  )

  heat_tree_pq(GPsubset,
    node_size = nb_sequences,
    node_color = n_obs,
    node_label = taxon_names,
    tree_label = taxon_names,
    node_size_trans = "log10 area"
  )
}
#> Loading required namespace: metacoder
#> This is metacoder version 0.3.7 (stable)
#> 
#> Attaching package: ‘metacoder’
#> The following object is masked from ‘package:ggplot2’:
#> 
#>     map_data
#> The following object is masked from ‘package:phyloseq’:
#> 
#>     filter_taxa
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’ ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’ ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’ ‘tidytree’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’ ‘tidytree’
#> Error in eval(e, x, parent.frame()): object 'GPsubset' not found
# }
```
