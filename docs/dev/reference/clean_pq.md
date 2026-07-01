# Clean phyloseq object by removing empty samples and taxa

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

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
  prefix_taxa_names = "_Taxa",
  check_taxonomy = FALSE,
  tax_remove_border_spaces = FALSE,
  tax_remove_all_space = FALSE,
  tax_replace_to_NA = FALSE,
  tax_redundant_suffix = FALSE,
  tax_replace_space_with = "_",
  tax_replace_invisible_chars = FALSE,
  tax_replace_NA_string = FALSE
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

  (logical, default FALSE) if TRUE, taxa (ASV, OTU) are renamed by their
  position in the OTU_table and prefix_taxa_names param (Taxa_1, Taxa_2,
  ...). Default to FALSE. If rename taxa (ASV, OTU) is true, the taxa
  (ASV, OTU) names in verbose information can be misleading.

- simplify_taxo:

  (logical) if TRUE, correct the taxonomy_table using the
  [`MiscMetabar::simplify_taxo()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/simplify_taxo.md)
  function

- prefix_taxa_names:

  (default "Taxa\_"): the prefix of taxa names (eg. "ASV\_" or "OTU\_")

- check_taxonomy:

  (logical, default FALSE) If TRUE, call
  [`verify_tax_table()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/verify_tax_table.md)
  to check for common taxonomy table issues.

- tax_remove_border_spaces:

  (logical, default FALSE) If TRUE, trim leading/trailing whitespace
  from values in the `tax_table` slot (passed to
  [`verify_tax_table()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/verify_tax_table.md)
  with `modify_phyloseq = TRUE`). Handles both ASCII whitespace and
  Unicode separators such as NBSP (U+00A0).

- tax_remove_all_space:

  (logical, default FALSE) If TRUE, replace internal whitespace (ASCII
  or Unicode separator) in `tax_table` values with `replace_space_with`.

- tax_replace_to_NA:

  (logical or character, default FALSE) If TRUE, replace `tax_table`
  values matching the default
  [unwanted_tax_patterns](https://adrientaudiere.github.io/MiscMetabar/dev/reference/unwanted_tax_patterns.md)
  with `NA`. A character vector of regex patterns can be supplied to
  override the defaults.

- tax_redundant_suffix:

  (logical or character, default FALSE) If TRUE, replace redundant
  `"_sp"` values with `NA` (e.g. `Russula_sp` at Species when `Russula`
  is already at Genus). A character string supplies a custom suffix.

- tax_replace_space_with:

  (character, default `"_"`) Replacement for internal whitespace when
  `tax_remove_all_space = TRUE`.

- tax_replace_invisible_chars:

  (logical, default FALSE) If TRUE, strip invisible / unusual characters
  (control chars, zero-width space, NBSP inside values, ...) from
  `tax_table` values. See
  [`verify_tax_table()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/verify_tax_table.md)'s
  `replace_invisible_chars` for the exact pattern.

- tax_replace_NA_string:

  (logical, default FALSE) If TRUE, replace the literal strings `"NA"`,
  `"NA NA"`, `"NA NA NA"` (any whitespace-separated repetition of `NA`,
  a common artifact of pasting taxonomic ranks together) in `tax_table`
  values with true `<NA>`. Case-sensitive to avoid clobbering real data.
  Default `FALSE` to avoid breaking changes.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object

## Author

Adrien Taudière

## Examples

``` r
clean_pq(data_fungi_mini)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]
# \donttest{
# Trim leading/trailing whitespace in tax_table values
clean_pq(data_fungi_mini, tax_remove_border_spaces = TRUE)
#> No values to modify. Returning original phyloseq object.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]
# Replace NA-like values (e.g. "unidentified", "NA") with NA
clean_pq(data_fungi_mini, tax_replace_to_NA = TRUE)
#> Replaced 4 NA-like value(s) with NA. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis
#> Total: 4 modification(s) in the taxonomy table: 4 value(s) replaced with NA (4 NA-like patterns, 0 short values, 0 redundant suffixes).
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]
# Drop redundant "_sp" tips
clean_pq(data_fungi_mini, tax_redundant_suffix = TRUE)
#> No values to modify. Returning original phyloseq object.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]
# Replace "NA" / "NA NA" concatenation artifacts with true <NA>
clean_pq(data_fungi_mini, tax_replace_NA_string = TRUE)
#> No values to modify. Returning original phyloseq object.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]
# }
```
