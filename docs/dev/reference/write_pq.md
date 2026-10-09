# Save phyloseq object in the form of multiple csv tables.

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This is the reverse function of
[`read_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/read_pq.md).

## Usage

``` r
write_pq(
  physeq,
  path = NULL,
  rdata = FALSE,
  one_file = FALSE,
  write_sam_data = TRUE,
  sam_data_first = FALSE,
  clean_pq = TRUE,
  reorder_taxa = FALSE,
  rename_taxa = FALSE,
  remove_empty_samples = TRUE,
  remove_empty_taxa = TRUE,
  clean_samples_names = TRUE,
  silent = FALSE,
  verbose = FALSE,
  quote = FALSE,
  sep_csv = "\t",
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- path:

  a path to the folder to save the phyloseq object

- rdata:

  (logical) does the phyloseq object is also saved in Rdata format?

- one_file:

  (logical) if TRUE, combine all data in one file only

- write_sam_data:

  (logical) does the samples data are add to the file. Only used if
  `one_file` is TRUE. Note that these option result in a lot of NA
  values.

- sam_data_first:

  (logical) if TRUE, put the sample data at the top of the table Only
  used if `one_file` and write_sam_data are both TRUE.

- clean_pq:

  (logical) If set to TRUE, empty samples are discarded after subsetting
  taxa (ASV, OTU, ...)

- reorder_taxa:

  (logical) if TRUE the otu_table is ordered by the number of sequences
  of taxa (ASV, OTU, ...) (descending order). Default to TRUE. Only
  possible if clean_pq is set to TRUE.

- rename_taxa:

  reorder_taxa (logical) if TRUE, taxa (ASV, OTU, ...) are renamed by
  their position in the OTU_table (asv_1, asv_2, ...). Default to FALSE.
  Only possible if clean_pq is set to TRUE.

- remove_empty_samples:

  (logical) Do you want to remove samples without sequences (this is
  done after removing empty taxa)

- remove_empty_taxa:

  (logical) Do you want to remove taxa without sequences (this is done
  before removing empty samples)

- clean_samples_names:

  (logical) Do you want to clean samples names?

- silent:

  (logical) If true, no message are printing.

- verbose:

  (logical) Additional informations in the message the verbose parameter
  overwrite the silent parameter.

- quote:

  a logical value (default FALSE) or a numeric vector. If TRUE, any
  character or factor columns will be surrounded by double quotes. If a
  numeric vector, its elements are taken as the indices of columns to
  quote. In both cases, row and column names are quoted if they are
  written. If FALSE nothing is quoted.

- sep_csv:

  (default tabulation) separator for column

- ...:

  Additional arguments passed on to
  [`utils::write.table()`](https://rdrr.io/r/utils/write.table.html)
  function.

## Value

Build a folder (path) containing one to four csv tables (refseq.csv,
otu_table.csv, tax_table.csv, sam_data.csv) and if present a phy_tree in
Newick format

## See also

[`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/save_pq.md)

## Author

Adrien Taudière

## Examples

``` r
write_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"))
write_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"), one_file = TRUE)
unlink(paste0(tempdir(), "/phyloseq"), recursive = TRUE)
```
