# Default patterns for unwanted taxonomic values

A named character vector of regular expressions used to identify common
problematic values in taxonomy tables. Each element is a regex pattern;
names provide human-readable descriptions.

Used as the default `replace_to_NA` argument in
[`verify_tax_table()`](https://adrientaudiere.github.io/MiscMetabar/reference/verify_tax_table.md)
and can be reused by other pqverse packages (e.g.
`dbpq::count_unwanted_tax()`).

## Usage

``` r
unwanted_tax_patterns
```

## Format

A named character vector with 17 elements:

- NA-like (NA, NaN, nan):

  `"^[Nn][Aa][Nn]?$"`

- NA-like (N/A, n/a):

  `"^[Nn]/[Aa]$"`

- None / none:

  `"^[Nn]one$"`

- empty string:

  `"^$"`

- whitespace only:

  `"^\\\\s+$"`

- unclassified:

  `"[Uu]nclassified"`

- unknown:

  `"[Uu]nknown"`

- unidentified:

  `"[Uu]nidentified"`

- uncultured:

  `"[Uu]ncultured"`

- incertae sedis:

  `"[Ii]ncertae[_\\\\s]?[Ss]edis"`

- metagenome:

  `"^[Mm]etagenome$"`

- environmental:

  `"^[Ee]nvironmental"`

- empty QIIME-style rank:

  `"^[kpcofgs]__$"`

- unknown species (\_sp prefix):

  `"^_sp"`

- unknown species (\_species prefix):

  `"^_species"`

- unknown cluster (MMseqs2):

  `"_uc$"`

- unknown ranks (PR2 database):

  `"__X+$"`

## See also

[`verify_tax_table()`](https://adrientaudiere.github.io/MiscMetabar/reference/verify_tax_table.md)

## Examples

``` r
unwanted_tax_patterns
#>                                 NA-like (NA, NaN, nan) 
#>                                      "^[Nn][Aa][Nn]?$" 
#>                                     NA-like (N/A, n/a) 
#>                                          "^[Nn]/[Aa]$" 
#>                                            None / none 
#>                                            "^[Nn]one$" 
#>                                           empty string 
#>                                                   "^$" 
#>                                        whitespace only 
#>                                               "^\\s+$" 
#>                                           unclassified 
#>                                      "[Uu]nclassified" 
#>                                                unknown 
#>                                           "[Uu]nknown" 
#>                                           unidentified 
#>                                      "[Uu]nidentified" 
#>                                             uncultured 
#>                                        "[Uu]ncultured" 
#>                                         incertae sedis 
#>                           "[Ii]ncertae[_\\s]?[Ss]edis" 
#>                                             metagenome 
#>                                      "^[Mm]etagenome$" 
#>                                          environmental 
#>                                    "^[Ee]nvironmental" 
#>                                 empty QIIME-style rank 
#>                                        "^[kpcofgs]__$" 
#>                           unknown species (_sp prefix) 
#>                                                 "^_sp" 
#>                      unknown species (_species prefix) 
#>                                            "^_species" 
#> unknown cluster (_uc prefix, e.g. MMseqs2 assignation) 
#>                                                 "_uc$" 
#>                           unknown ranks (PR2 database) 
#>                                                "__X+$" 
# Use with grepl to check a value
any(vapply(
  unwanted_tax_patterns,
  \(pat) grepl(pat, "unclassified"),
  logical(1)
))
#> [1] TRUE
```
