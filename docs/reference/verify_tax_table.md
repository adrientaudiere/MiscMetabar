# Verify the taxonomy table of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Check taxonomy table for common issues and send warnings/messages
accordingly. This function is called by
[`verify_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/verify_pq.md)
when `check_taxonomy = TRUE`.

## Usage

``` r
verify_tax_table(
  physeq,
  verbose = FALSE,
  replace_to_NA = c("^[Nn][Aa][Nn]?$", "^[Nn]/[Aa]$", "^[Nn]one$", "^$", "^\\s+$",
    "[Uu]nclassified", "[Uu]nknown", "[Uu]nidentified", "[Uu]ncultured",
    "[Ii]ncertae[_\\s]?[Ss]edis", "^[Mm]etagenome$", "^[Ee]nvironmental",
    "^[kpcofgs]__$"),
  min_char = 4,
  redundant_suffix = "_sp",
  taxonomic_ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  modify_phyloseq = FALSE,
  remove_border_spaces = TRUE,
  remove_all_space = FALSE,
  replace_space_with = "_"
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- verbose:

  (logical, default FALSE) If TRUE, print warnings and messages about
  potential taxonomy issues.

- replace_to_NA:

  (character vector) A vector of regex patterns to identify values that
  should be considered as NA. Default patterns include common
  placeholders like "unclassified", "unknown", "uncultured",
  "incertae_sedis", "metagenome", empty QIIME-style ranks (e.g.,
  "k\_\_"), etc.

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

  - Leading/trailing whitespace (if `remove_border_spaces = TRUE`)

  - Internal spaces (if `remove_all_space = TRUE`)

  Messages will indicate the number of values replaced for each type.

- remove_border_spaces:

  (logical, default TRUE) If TRUE and `modify_phyloseq = TRUE`, remove
  leading and trailing whitespace from taxonomic values.

- remove_all_space:

  (logical, default FALSE) If TRUE and `modify_phyloseq = TRUE`, replace
  internal spaces (spaces within taxonomic values) with the character
  specified in `replace_space_with`.

- replace_space_with:

  (character, default "\_") Character to use when replacing internal
  spaces. Only used when `remove_all_space = TRUE`.

## Value

If `modify_phyloseq = FALSE` (default): Nothing (invisible NULL).
Warnings/messages only if verbose = TRUE and issues are found. If
`modify_phyloseq = TRUE`: The modified phyloseq object with problematic
values replaced by NA, along with messages summarizing the changes.

## Author

Adrien Taudière

## Examples

``` r
verify_tax_table(data_fungi)
# \donttest{
verify_tax_table(data_fungi, verbose = TRUE)
#> Warning: Found 145 taxonomic value(s) matching NA-like patterns. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, .... Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 1432 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 976 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1319 taxonomic value(s) with internal spaces: 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Check for redundant "_sp" patterns (default)
data_fungi2 <- data_fungi
data_fungi2@tax_table[1, "Species"] <- "Eutypa_sp"
verify_tax_table(data_fungi2, verbose = TRUE, redundant_suffix = "_sp")
#> Warning: Found 145 taxonomic value(s) matching NA-like patterns. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, .... Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 1432 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 975 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1319 taxonomic value(s) with internal spaces: 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.
#> Warning: Found 1 taxonomic value(s) with redundant '_sp' patterns where the information is already present at a higher rank: Eutypa_sp (Species = Genus). Use modify_phyloseq = TRUE to replace these with NA.

# Automatically replace problematic values with NA
# This replaces: NA-like patterns, short values, and redundant suffixes
data_fungi2_cleaned <- verify_tax_table(data_fungi2,
  modify_phyloseq = TRUE
)
#> Replaced 145 NA-like value(s) with NA. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, ...
#> Replaced 1432 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Replaced 1 redundant '_sp' value(s) with NA: Eutypa_sp (Species = Genus)
#> Total: 1578 modification(s) in the taxonomy table: 1578 value(s) replaced with NA (145 NA-like patterns, 1432 short values, 1 redundant suffixes).
# Check that the redundant value was replaced
data_fungi2@tax_table[1, "Species"] # "Eutypa_sp"
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Species    
#> ASV2 "Eutypa_sp"
data_fungi2_cleaned@tax_table[1, "Species"] # NA
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Species
#> ASV2 NA     

# Combine verbose mode with modifications to see all issues
data_fungi2_cleaned <- verify_tax_table(data_fungi2,
  verbose = TRUE,
  modify_phyloseq = TRUE
)
#> Replaced 145 NA-like value(s) with NA. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, ...
#> Replaced 1432 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Found 976 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1319 taxonomic value(s) with internal spaces: 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.
#> Replaced 1 redundant '_sp' value(s) with NA: Eutypa_sp (Species = Genus)
#> Total: 1578 modification(s) in the taxonomy table: 1578 value(s) replaced with NA (145 NA-like patterns, 1432 short values, 1 redundant suffixes).

# Check for other patterns like "_var" or "_cf"
verify_tax_table(data_fungi, verbose = TRUE, redundant_suffix = "_var")
#> Warning: Found 145 taxonomic value(s) matching NA-like patterns. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, .... Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 1432 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 976 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1319 taxonomic value(s) with internal spaces: 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Disable redundant suffix check
verify_tax_table(data_fungi, verbose = TRUE, redundant_suffix = NULL)
#> Warning: Found 145 taxonomic value(s) matching NA-like patterns. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, .... Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 1432 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 976 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1319 taxonomic value(s) with internal spaces: 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Specify custom taxonomic rank order
verify_tax_table(data_fungi,
  verbose = TRUE,
  taxonomic_ranks = c("Class", "Order", "Family", "Genus")
)
#> Warning: Found 145 taxonomic value(s) matching NA-like patterns. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, .... Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 1432 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 976 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1319 taxonomic value(s) with internal spaces: 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Handle whitespace in taxonomic values
# Create example with spaces
data_fungi3 <- data_fungi
data_fungi3@tax_table[1, "Genus"] <- " Russula "
data_fungi3@tax_table[2, "Species"] <- "Russula emetica"

# Check for spaces (verbose mode)
verify_tax_table(data_fungi3, verbose = TRUE)
#> Warning: Found 145 taxonomic value(s) matching NA-like patterns. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, .... Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 1432 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 974 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 1 taxonomic value(s) with leading or trailing whitespace: ' Russula ' (Genus). Use modify_phyloseq = TRUE to trim these values.
#> Warning: Found 1320 taxonomic value(s) with internal spaces: 'Russula emetica' (Species), 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.

# Remove leading/trailing whitespace (enabled by default)
data_fungi3_trimmed <- verify_tax_table(data_fungi3, modify_phyloseq = TRUE)
#> Replaced 145 NA-like value(s) with NA. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, ...
#> Replaced 1432 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Trimmed leading/trailing whitespace from 1 value(s): ' Russula ' (Genus)
#> Total: 1578 modification(s) in the taxonomy table: 1577 value(s) replaced with NA (145 NA-like patterns, 1432 short values, 0 redundant suffixes); 1 value(s) trimmed of border whitespace.
data_fungi3_trimmed@tax_table[1, "Genus"] # "Russula" (trimmed)
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Genus    
#> ASV2 "Russula"

# Also replace internal spaces with underscores
data_fungi3_cleaned <- verify_tax_table(data_fungi3,
  modify_phyloseq = TRUE,
  remove_all_space = TRUE,
  replace_space_with = "_"
)
#> Replaced 145 NA-like value(s) with NA. Unique values: Pezizomycotina_cls_Incertae_sedis, Rozellomycotina_cls_Incertae_sedis, Sordariomycetes_ord_Incertae_sedis, Microbotryomycetes_ord_Incertae_sedis, Lecanoromycetes_ord_Incertae_sedis, Dothideomycetes_ord_Incertae_sedis, Cystobasidiomycetes_ord_Incertae_sedis, Pezizomycotina_ord_Incertae_sedis, Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, ...
#> Replaced 1432 short value(s) (< 4 chars) with NA: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars)
#> Trimmed leading/trailing whitespace from 1 value(s): ' Russula ' (Genus)
#> Replaced internal spaces with '_' in 1320 value(s): 'Russula emetica' (Species), 'Plant Pathogen' (Guild), 'Endophyte-Undefined Saprotroph-Wood Saprotroph' (Guild), 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), ...
#> Total: 2898 modification(s) in the taxonomy table: 1577 value(s) replaced with NA (145 NA-like patterns, 1432 short values, 0 redundant suffixes); 1 value(s) trimmed of border whitespace; 1320 value(s) had internal spaces replaced.
data_fungi3_cleaned@tax_table[2, "Species"] # "Russula_emetica"
#> Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
#>      Species          
#> ASV6 "Russula_emetica"
# }
```
