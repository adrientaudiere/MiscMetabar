# Pseudo-merge taxa in groups

Returns `x` pruned to the first taxon of each group defined in `group`.

## Usage

``` r
merge_taxa_vec_pseudo(x, group, reorder = FALSE)
```

## Arguments

- x:

  a phyloseq component-class object

- group:

  a vector with one element for each taxon in `x` that defines the new
  groups
