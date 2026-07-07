# Export a phyloseq object to GBIF MDT template CSV/TSV files

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Write the OTU table, the taxonomy table and the sample data of a
phyloseq object to separate delimited text files, one per table,
following the GBIF [Metabarcoding Data Toolkit
(MDT)](https://mdt.gbif.org/) template layout. This is the plain-text
(CSV/TSV) counterpart of
[`phyloseq_to_MDT_excel()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/phyloseq_to_MDT_excel.md),
which writes a single multi-sheet `.xlsx` file.

The MDT template expects the OTU table with **OTU IDs in rows and sample
IDs in columns** (sequence read counts in the cells), and the sample and
taxonomy tables keyed by a leading `id` column. The reference sequences,
when available, are appended as a `DNA_sequence` column of the taxonomy
file (Darwin Core DNA-derived-data extension term).

When `check_dwc = TRUE` (the default), a lightweight Darwin Core
compliance check warns about recommended sample-level terms that are
missing from `sample_data` (`decimalLatitude`, `decimalLongitude`,
`eventDate`). This is a non-blocking helper, not a full validation of
the GBIF MDT template.

## Usage

``` r
phyloseq_to_MDT_csv(
  physeq,
  path = ".",
  prefix = "",
  sep = ",",
  check_dwc = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object.

- path:

  (character, default `"."`) Directory where the files are written.
  Created (recursively) if it does not exist.

- prefix:

  (character, default `""`) Optional prefix prepended to each output
  file name (e.g. `"data_fungi_"`).

- sep:

  (character, default `","`) Field separator. Use `"\t"` to write
  tab-separated (TSV) files, the format favored by the GBIF MDT
  validator. The file extension follows `sep` (`.csv` for `","`, `.tsv`
  otherwise).

- check_dwc:

  (logical, default TRUE) If TRUE, warn about recommended Darwin Core
  sample-level terms missing from `sample_data`.

## Value

Invisibly returns a named character vector of the written file paths
(`OTU_table`, `Taxonomy`, `Samples`).

## Details

See the GBIF [Metabarcoding Data Toolkit](https://mdt.gbif.org/) for the
expected input format. The MDT accepts both TSV and XLSX uploads.

## See also

[`phyloseq_to_MDT_excel()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/phyloseq_to_MDT_excel.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
out_dir <- file.path(tempdir(), "mdt_csv")
files <- phyloseq_to_MDT_csv(clean_pq(data_fungi_mini), path = out_dir)
#> Warning: ! Recommended Darwin Core sample terms missing from sample_data:
#>   "decimalLatitude", "decimalLongitude", and "eventDate".
#> ℹ Add them before GBIF MDT submission if available.
#> MDT template files written to /tmp/RtmprMwnI8/mdt_csv
file.exists(files)
#> [1] TRUE TRUE TRUE
unlink(out_dir, recursive = TRUE)
# }
if (FALSE) { # \dontrun{
# Tab-separated output (MDT-favored TSV)
phyloseq_to_MDT_csv(data_fungi_mini, path = "mdt", sep = "\t")
} # }
```
