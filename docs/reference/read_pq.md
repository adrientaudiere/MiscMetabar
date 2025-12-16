# Read phyloseq object from multiple csv tables and a phylogenetic tree in Newick format.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This is the reverse function of
[`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/write_pq.md).

## Usage

``` r
read_pq(
  path = NULL,
  taxa_are_rows = FALSE,
  sam_names = NULL,
  sep_csv = "\t",
  ...
)
```

## Arguments

- path:

  (required) a path to the folder to read the phyloseq object

- taxa_are_rows:

  (default to FALSE) see ?phyloseq for details

- sam_names:

  The name of the variable (column) in sam_data.csv to rename samples.
  Note that if you use
  [`write_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/MiscMetabar-deprecated.md)
  function to save your physeq object, you may use sam_names = "X" to
  rename the samples names as before.

- sep_csv:

  (default tabulation) separator for column

- ...:

  Additional arguments passed on to
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html)
  function.

## Value

One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv,
sam_data.csv) and if present a phy_tree in Newick format. At least the
otu_table.csv need to be present.

## Author

Adrien Taudière

## Examples

``` r
write_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"))
read_pq(path = paste0(tempdir(), "/phyloseq"))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 1420 taxa and 185 samples ]
#> sample_data() Sample Data:       [ 185 samples by 8 sample variables ]
#> tax_table()   Taxonomy Table:    [ 1420 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 1420 reference sequences ]
unlink(paste0(tempdir(), "/phyloseq"), recursive = TRUE)
```
