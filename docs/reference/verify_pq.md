# Verify the validity of a phyloseq object

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mostly for internal use in MiscMetabar functions.

## Usage

``` r
verify_pq(
  physeq,
  verbose = FALSE,
  min_nb_seq_sample = 500,
  min_nb_seq_taxa = 1,
  check_taxonomy = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- verbose:

  (logical, default FALSE) If TRUE, prompt some warnings.

- min_nb_seq_sample:

  (numeric) Only used if verbose = TRUE. Minimum number of sequences per
  samples to not show warning.

- min_nb_seq_taxa:

  (numeric) Only used if verbose = TRUE. Minimum number of sequences per
  taxa to not show warning.

- check_taxonomy:

  (logical, default FALSE) If TRUE, call
  [`verify_tax_table()`](https://adrientaudiere.github.io/MiscMetabar/reference/verify_tax_table.md)
  to check for common taxonomy table issues.

- ...:

  Additional arguments passed to
  [`verify_tax_table()`](https://adrientaudiere.github.io/MiscMetabar/reference/verify_tax_table.md)
  when `check_taxonomy = TRUE`.

## Value

Nothing if the phyloseq object is valid. An error in the other case.
Warnings if verbose = TRUE or check_taxonomy = TRUE

## Author

Adrien Taudière

## Examples

``` r

verify_pq(data_fungi_mini)
# \donttest{
verify_pq(data_fungi_mini, check_taxonomy = TRUE)
#> Warning: Found 4 taxonomic value(s) matching NA-like patterns. Unique values: Cantharellales_fam_Incertae_sedis, Atractiellales_fam_Incertae_sedis, Russulales_fam_Incertae_sedis, Hymenochaetales_fam_Incertae_sedis. Use modify_phyloseq = TRUE to replace these with NA.
#> Warning: Found 20 taxonomic value(s) with less than 4 characters: - (Trophic.Mode, 1 chars), - (Guild, 1 chars), - (Trait, 1 chars), - (Confidence.Ranking, 1 chars). Use modify_phyloseq = TRUE to replace these with NA.
#> Found 13 taxa with duplicate taxonomic paths. This may indicate redundant taxa or issues with taxonomic assignment.
#> Warning: Found 70 taxonomic value(s) with internal spaces: 'Wood Saprotroph-Undefined Saprotroph' (Guild), 'Undefined Saprotroph' (Guild), 'Wood Saprotroph' (Guild), 'Ectomycorrhizal-Wood Saprotroph' (Guild), 'Leaf Saprotroph-Plant Pathogen-Undefined Saprotroph-Wood Saprotroph' (Guild), .... Use modify_phyloseq = TRUE and remove_all_space = TRUE to replace these spaces with '_'.
# }
```
