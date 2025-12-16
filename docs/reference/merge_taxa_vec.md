# Merge taxa in groups (vectorized version)

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Firstly release in the [speedyseq](https://github.com/mikemc/speedyseq/)
R package by Michael R. McLaren.

Merge taxa in `x` into a smaller set of taxa defined by the vector
`group`. Taxa whose value in `group` is NA will be dropped. New taxa
will be named according to the most abundant taxon in each group
(`phyloseq` and `otu_table` objects) or the first taxon in each group
(all other phyloseq component objects).

If `x` is a phyloseq object with a phylogenetic tree, then the new taxa
will be ordered as they are in the tree. Otherwise, the taxa order can
be controlled by the `reorder` argument, which behaves like the
`reorder` argument in
[`base::rowsum()`](https://rdrr.io/r/base/rowsum.html).
`reorder = FALSE` will keep taxa in the original order determined by
when the member of each group first appears in `taxa_names(x)`;
`reorder = TRUE` will order new taxa according to their corresponding
value in `group`.

The `tax_adjust` argument controls the handling of taxonomic
disagreements within groups. Setting `tax_adjust == 0` causes no
adjustment; the taxonomy of the new group is set to the archetype taxon
(the most abundant taxon in each group). Otherwise, disagreements within
a group at a given rank cause the values at lower ranks to be set to
`NA`. If `tax_adjust == 1` (the default), then a rank where all taxa in
the group are already NA is not counted as a disagreement, and lower
ranks may be kept if the taxa agree. This corresponds to the original
phyloseq behavior. If `tax_adjust == 2`, then these NAs are treated as a
disagreement; all ranks are set to NA after the first disagreement or
NA.

## Usage

``` r
merge_taxa_vec(x, group, reorder = FALSE, tax_adjust = 1L)

# S4 method for class 'phyloseq'
merge_taxa_vec(x, group, reorder = FALSE, tax_adjust = 1L)

# S4 method for class 'otu_table'
merge_taxa_vec(x, group, reorder = FALSE)

# S4 method for class 'taxonomyTable'
merge_taxa_vec(x, group, reorder = FALSE, tax_adjust = 1L)

# S4 method for class 'phylo'
merge_taxa_vec(x, group)

# S4 method for class 'XStringSet'
merge_taxa_vec(x, group, reorder = FALSE)
```

## Arguments

- x:

  A phyloseq object or component object

- group:

  A vector with one element for each taxon in `physeq` that defines the
  new groups. see
  [`base::rowsum()`](https://rdrr.io/r/base/rowsum.html).

- reorder:

  Logical specifying whether to reorder the taxa by their `group`
  values. Ignored if `x` has (or is) a phylogenetic tree.

- tax_adjust:

  0: no adjustment; 1: phyloseq-compatible adjustment; 2: conservative
  adjustment

## Value

A new phyloseq-class, otu_table, tax_table, XStringset or sam_data
object depending on the class of the x param

## See also

Function in MiscMetabar that use this function:
[`postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md)

[`base::rowsum()`](https://rdrr.io/r/base/rowsum.html)

[`phyloseq::merge_taxa()`](https://rdrr.io/pkg/phyloseq/man/merge_taxa-methods.html)

## Author

Michael R. McLaren (orcid:
[0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X)) modified
by Adrien Taudiere
