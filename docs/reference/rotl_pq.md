# rotl wrapper for phyloseq data

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Make a taxonomic tree using the ASV names of a physeq object and the
Open Tree of Life tree.

## Usage

``` r
rotl_pq(
  physeq,
  taxonomic_rank = c("Genus", "Species"),
  context_name = "All life",
  discard_genus_alone = TRUE,
  pattern_to_remove_tip = c("ott\\d+|_ott\\d+"),
  pattern_to_remove_node = c("_ott.*|mrca*")
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxonomic_rank:

  (Character) The column(s) present in the @tax_table slot of the
  phyloseq object. Can be a vector of two columns (e.g. the default
  c("Genus", "Species")). If only one column is set it need to be format
  in this way ("Genus species" for ex. "Quercus robur") with a space.

- context_name:

  : can bue used to select only a part of the Open Tree of Life. See
  `?rotl::tnrs_contexts()` for available values

- discard_genus_alone:

  (logical) If TRUE (default), genus without information at the species
  level are discarded.

- pattern_to_remove_tip:

  (character regex string) A regex to remove unwanted part of tip names.
  If set to null, tip names are left intact.

- pattern_to_remove_node:

  (character regex string) A regex to remove unwanted part of node
  names. If set to null, node names are left intact.

## Value

A plot

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to `rotl` package if you use this function.

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("rotl")) {
  tr <- rotl_pq(data_fungi_mini, pattern_to_remove_tip = NULL)
  plot(tr)

  tr_Asco <- rotl_pq(data_fungi,
    taxonomic_rank = c("Genus", "Species"),
    context_name = "Ascomycetes"
  )
  plot(tr_Asco)
}
} # }
```
