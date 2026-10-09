# A wrapper of write_pq to save in all three possible formats

A wrapper of write_pq to save in all three possible formats

## Usage

``` r
save_pq(physeq, path = NULL, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- path:

  a path to the folder to save the phyloseq object

- ...:

  Additional arguments passed on to
  [`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/write_pq.md)
  or [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html)
  function.

## Value

Build a folder (in path) with four csv tables (`refseq.csv`,
`otu_table.csv`, `tax_table.csv`, `sam_data.csv`) + one table with all
tables together + a rdata file (`physeq.RData`) that can be loaded using
[`base::load()`](https://rdrr.io/r/base/load.html) function + if present
a phylogenetic tree in Newick format (`phy_tree.txt`)

## Details

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Write :

- 4 separate tables

- 1 table version

- 1 RData file

## See also

[`write_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/write_pq.md)

## Author

Adrien Taudière

## Examples

``` r
save_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"))
unlink(paste0(tempdir(), "/phyloseq"), recursive = TRUE)
```
