# Export a phyloseq object to a multi-sheet Excel file for GBIF MDT submission

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Write the OTU table, the sample data, the taxonomy table and (if
present) the reference sequences of a phyloseq object to a single
multi-sheet `.xlsx` file, formatted for submission to the GBIF
[Metabarcoding Data Toolkit (MDT)](https://mdt.gbif.org/).

Each sheet gets a leading `id` column holding the row identifiers (taxa
or samples names). The reference sequences, when available, are appended
as a `DNA_sequence` column of the taxonomy sheet (Darwin Core
DNA-derived-data extension term).

When `check_dwc = TRUE` (the default), a lightweight Darwin Core
compliance check warns about recommended sample-level terms that are
missing from `sample_data` (`decimalLatitude`, `decimalLongitude`,
`eventDate`). This is a non-blocking helper, not a full validation of
the GBIF MDT template.

## Usage

``` r
phyloseq_to_MDT_excel(
  physeq,
  filename = "Phyloseq_Tables.xlsx",
  check_dwc = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object.

- filename:

  (character, default `"Phyloseq_Tables.xlsx"`) Path of the output
  `.xlsx` file.

- check_dwc:

  (logical, default TRUE) If TRUE, warn about recommended Darwin Core
  sample-level terms missing from `sample_data`.

## Value

Invisibly returns the path to the written file (`filename`).

## Details

This function requires the writexl package. See the GBIF [Metabarcoding
Data Toolkit](https://mdt.gbif.org/) for the expected input format.

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("writexl")) {
  out <- file.path(tempdir(), "data_fungi_mini_MDT.xlsx")
  phyloseq_to_MDT_excel(clean_pq(data_fungi_mini), filename = out)
  file.exists(out)
  unlink(out)
}
#> Loading required namespace: writexl
#> Warning: ! Recommended Darwin Core sample terms missing from sample_data:
#>   "decimalLatitude", "decimalLongitude", and "eventDate".
#> ℹ Add them before GBIF MDT submission if available.
#> Excel file written to /tmp/RtmplwIHEj/data_fungi_mini_MDT.xlsx
# }
```
