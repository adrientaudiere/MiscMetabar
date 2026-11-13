# Verify the taxonomy table of a phyloseq object

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Check taxonomy table for common issues and send warnings/messages
accordingly. This function is called by
[`verify_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/verify_pq.md)
when `check_taxonomy = TRUE`.

## Usage

``` r
verify_tax_table(
  physeq,
  verbose = TRUE,
  replace_to_NA = unwanted_tax_patterns,
  min_char = 4,
  redundant_suffix = "_sp",
  taxonomic_ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  modify_phyloseq = FALSE,
  remove_border_spaces = TRUE,
  remove_all_space = FALSE,
  replace_space_with = "_",
  detect_invisible_chars = TRUE,
  replace_invisible_chars = FALSE,
  invisible_chars_replacement = ""
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- verbose:

  (logical, default TRUE) If TRUE, print warnings and messages about
  potential taxonomy issues.

- replace_to_NA:

  (character vector) A vector of regex patterns to identify values that
  should be considered as NA. Defaults to
  [unwanted_tax_patterns](https://adrientaudiere.github.io/MiscMetabar/dev/reference/unwanted_tax_patterns.md),
  a named character vector of common placeholders like "unclassified",
  "unknown", "uncultured", "incertae_sedis", "metagenome", empty
  QIIME-style ranks, etc.

- min_char:

  (integer, default 4) Minimum number of characters for a taxonomic
  value to be considered valid. Values with fewer characters (excluding
  NA) will trigger a warning when verbose = TRUE.

- redundant_suffix:

  (character, default "\_sp") Suffix pattern to detect redundant
  taxonomic information. For example, "Russula_sp" in Species column is
  redundant if "Russula" is already present in the Genus column. Set to
  NULL to disable this check. Other examples: "\_var", "\_ssp", "\_cf".

- taxonomic_ranks:

  (character vector, default NULL) Names of taxonomic ranks in
  hierarchical order from highest to lowest (e.g., c("Kingdom",
  "Phylum", "Class", "Order", "Family", "Genus", "Species")). If NULL,
  uses the column names of the taxonomy table in their existing order.
  Used to determine parent-child relationships for redundant suffix
  detection.

- modify_phyloseq:

  (logical, default FALSE) If TRUE, replace problematic values with NA
  in the taxonomy table and return the modified phyloseq object. The
  following types of values are replaced:

  - Values matching `replace_to_NA` patterns (e.g., "unclassified",
    "unknown")

  - Values with fewer than `min_char` characters

  - Redundant suffix patterns (e.g., "Russula_sp" when "Russula" is in
    Genus)

  - Leading/trailing whitespace, including non-breaking space U+00A0 and
    other Unicode separators (if `remove_border_spaces = TRUE`)

  - Internal spaces, including Unicode separators (if
    `remove_all_space = TRUE`)

  - Invisible / unusual characters such as control chars, zero-width
    space U+200B or non-breaking space U+00A0 inside values (if
    `replace_invisible_chars = TRUE`)

  Messages will indicate the number of values replaced for each type.

- remove_border_spaces:

  (logical, default TRUE) If TRUE and `modify_phyloseq = TRUE`, remove
  leading and trailing whitespace from taxonomic values. Matches both
  ASCII whitespace and Unicode separators (NBSP, em space, ideographic
  space, ...) — [`trimws()`](https://rdrr.io/r/base/trimws.html) alone
  only handles `[ \\t\\r\\n]` and would silently leave NBSP in place.

- remove_all_space:

  (logical, default FALSE) If TRUE and `modify_phyloseq = TRUE`, replace
  internal whitespace (ASCII or Unicode separator) with the character
  specified in `replace_space_with`.

- replace_space_with:

  (character, default "\_") Character to use when replacing internal
  spaces. Only used when `remove_all_space = TRUE`.

- detect_invisible_chars:

  (logical, default TRUE) If TRUE, scan taxonomic values for invisible /
  unusual characters: anything in Unicode category `\\p{C}` (control /
  format / surrogate / private use / unassigned) or any `\\p{Z}`
  separator other than a plain ASCII space or tab. Typical offenders
  include non-breaking space (U+00A0), zero-width space (U+200B),
  zero-width joiner (U+200D), and control characters. Letters with
  diacritics, digits and punctuation are NOT flagged.

- replace_invisible_chars:

  (logical, default FALSE) If TRUE and `modify_phyloseq = TRUE`, strip
  the characters detected by `detect_invisible_chars` from taxonomic
  values (replacement is `invisible_chars_replacement`, default empty
  string). Values that become empty after stripping are turned into
  `NA`.

- invisible_chars_replacement:

  (character, default `""`) Replacement string for
  `replace_invisible_chars = TRUE`.

## Value

If `modify_phyloseq = FALSE` (default): Nothing (invisible NULL).
Warnings/messages only if verbose = TRUE and issues are found. If
`modify_phyloseq = TRUE`: The modified phyloseq object with problematic
values replaced by NA, along with messages summarizing the changes.

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
verify_tax_table(data_fungi_mini)
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 13 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.
# }
```
