# Clean phyloseq object by removing empty samples and taxa

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

In addition, this function check for discrepancy (and rename) between
(i) taxa names in refseq, taxonomy table and otu_table and between (ii)
sample names in sam_data and otu_table.

## Usage

``` r
clean_pq(
  physeq,
  remove_empty_samples = TRUE,
  remove_empty_taxa = TRUE,
  clean_samples_names = TRUE,
  silent = FALSE,
  verbose = FALSE,
  force_taxa_as_columns = FALSE,
  force_taxa_as_rows = FALSE,
  reorder_taxa = FALSE,
  rename_taxa = FALSE,
  simplify_taxo = FALSE,
  prefix_taxa_names = "_Taxa"
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

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

- force_taxa_as_columns:

  (logical) If true, if the taxa are rows transpose the otu_table and
  set taxa_are_rows to false

- force_taxa_as_rows:

  (logical) If true, if the taxa are columns transpose the otu_table and
  set taxa_are_rows to true

- reorder_taxa:

  (logical) if TRUE the otu_table is ordered by the number of sequences
  of taxa (ASV, OTU) in descending order. Default to FALSE.

- rename_taxa:

  (logical) if TRUE, taxa (ASV, OTU) are renamed by their position in
  the OTU_table and prefix_taxa_names param (by default: Taxa_1, Taxa_2,
  ...). Default to FALSE. If rename taxa (ASV, OTU) is true, the taxa
  (ASV, OTU) names in verbose information can be misleading.

- simplify_taxo:

  (logical) if TRUE, correct the taxonomy_table using the
  [`MiscMetabar::simplify_taxo()`](https://adrientaudiere.github.io/MiscMetabar/reference/simplify_taxo.md)
  function

- prefix_taxa_names:

  (default "Taxa\_"): the prefix of taxa names (eg. "ASV\_" or "OTU\_")

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object

## Author

Adrien Taudière
