# Make Krona files using [KronaTools](https://github.com/marbl/Krona/wiki).

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Need the installation of kronatools on the computer ([installation
instruction](https://github.com/marbl/Krona/wiki/Installing)).

## Usage

``` r
krona(
  physeq,
  file_path = "krona.html",
  nb_seq = TRUE,
  ranks = "All",
  add_unassigned_rank = 0,
  name = NULL
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- file_path:

  (required) the location of the html file to save

- nb_seq:

  (logical) If true, Krona set the distribution of sequences in the
  taxonomy. If False, Krona set the distribution of ASVs in the
  taxonomy.

- ranks:

  Number of the taxonomic ranks to plot (num of the column in
  `tax_table` slot of your `physeq` object). Default setting plot all
  the ranks (argument 'All').

- add_unassigned_rank:

  (int) Add unassigned for rank inferior to 'add_unassigned_rank' when
  necessary.

- name:

  A name for intermediary files, Useful to name your krona result files
  before merging using
  [`merge_krona()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_krona.md).
  Must not contain space.

## Value

A html file

## Details

This function is mainly a wrapper of the work of others. Please cite
[Krona](https://github.com/marbl/Krona) if you use this function.

## See also

[`merge_krona`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_krona.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # tolower(Sys.info()[["sysname"]]) != "windows" && MiscMetabar::is_krona_installed()
data("GlobalPatterns", package = "phyloseq")
GA <- subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
if (FALSE) { # \dontrun{
krona(GA, "Number.of.sequences.html")
krona(GA, "Number.of.ASVs.html", nb_seq = FALSE)
merge_krona(c("Number.of.sequences.html", "Number.of.ASVs.html"))
} # }
}
```
