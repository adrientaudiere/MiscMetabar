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
verify_tax_table(data_fungi_mini, verbose = TRUE)
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 13 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Check for redundant "_sp" patterns (default)
data_fungi2 <- data_fungi_mini
data_fungi2@tax_table[1, "Species"] <- "Eutypa_sp"
verify_tax_table(data_fungi2, verbose = TRUE, redundant_suffix = "_sp")
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 12 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Automatically replace problematic values with NA
# This replaces: NA-like patterns, short values, and redundant suffixes
data_fungi2_cleaned <- verify_tax_table(data_fungi2,
  modify_phyloseq = TRUE
)
#> Replaced 4 NA-like value(s) with NA. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis
#> Replaced 20 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Found 12 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.
#> Total: 24 modification(s) in the taxonomy table: 24 value(s) replaced with NA (4 NA-like patterns, 20 short values, 0 redundant suffixes).
# Check that the redundant value was replaced
data_fungi2@tax_table[1, "Species"] # "Eutypa_sp"
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Species    
#> ASV7 "Eutypa_sp"
data_fungi2_cleaned@tax_table[1, "Species"] # NA
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Species    
#> ASV7 "Eutypa_sp"

# Combine verbose mode with modifications to see all issues
data_fungi2_cleaned <- verify_tax_table(data_fungi2,
  verbose = TRUE,
  modify_phyloseq = TRUE
)
#> Replaced 4 NA-like value(s) with NA. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis
#> Replaced 20 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Found 12 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.
#> Total: 24 modification(s) in the taxonomy table: 24 value(s) replaced with NA (4 NA-like patterns, 20 short values, 0 redundant suffixes).

# Check for other patterns like "_var" or "_cf"
verify_tax_table(data_fungi_mini, verbose = TRUE, redundant_suffix = "_var")
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 13 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Disable redundant suffix check
verify_tax_table(data_fungi_mini, verbose = TRUE, redundant_suffix = NULL)
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 13 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Specify custom taxonomic rank order
verify_tax_table(data_fungi_mini,
  verbose = TRUE,
  taxonomic_ranks = c("Class", "Order", "Family", "Genus")
)
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 13 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Handle whitespace in taxonomic values
# Create example with spaces
data_fungi3 <- data_fungi_mini
data_fungi3@tax_table[1, "Genus"] <- " Russula "
data_fungi3@tax_table[2, "Species"] <- "Russula emetica"

# Check for spaces (verbose mode)
verify_tax_table(data_fungi3, verbose = TRUE)
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 11 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1 taxonomic value(s) with leading or trailing whitespace: ' Russula ' (Genus). Use modify_phyloseq = TRUE to trim these values.
#> Warning: Found 71 taxonomic value(s) with internal spaces: 'Russula emetica' (Species), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Remove leading/trailing whitespace (enabled by default)
data_fungi3_trimmed <- verify_tax_table(data_fungi3, modify_phyloseq = TRUE)
#> Replaced 4 NA-like value(s) with NA. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis
#> Replaced 20 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Found 11 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Trimmed leading/trailing whitespace from 1 value(s): ' Russula ' (Genus)
#> Warning: Found 71 taxonomic value(s) with internal spaces: 'Russula emetica' (Species), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.
#> Total: 25 modification(s) in the taxonomy table: 24 value(s) replaced with NA (4 NA-like patterns, 20 short values, 0 redundant suffixes); 1 value(s) trimmed of border whitespace.
data_fungi3_trimmed@tax_table[1, "Genus"] # "Russula" (trimmed)
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Genus    
#> ASV7 "Russula"

# Also replace internal spaces with underscores
data_fungi3_cleaned <- verify_tax_table(data_fungi3,
  modify_phyloseq = TRUE,
  remove_all_space = TRUE,
  replace_space_with = "_"
)
#> Replaced 4 NA-like value(s) with NA. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis
#> Replaced 20 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Found 11 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Trimmed leading/trailing whitespace from 1 value(s): ' Russula ' (Genus)
#> Replaced internal spaces with '_' in 71 value(s): 'Russula emetica' (Species), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), ...
#> Total: 96 modification(s) in the taxonomy table: 24 value(s) replaced with NA (4 NA-like patterns, 20 short values, 0 redundant suffixes); 1 value(s) trimmed of border whitespace; 71 value(s) had internal spaces replaced.
data_fungi3_cleaned@tax_table[2, "Species"] # "Russula_emetica"
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Species          
#> ASV8 "Russula_emetica"
# }
```
