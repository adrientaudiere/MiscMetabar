# Merge Krona files using [KronaTools](https://github.com/marbl/Krona/wiki).

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Need the installation of kronatools on the computer ([installation
instruction](https://github.com/marbl/Krona/wiki/Installing)).

Function merge_krona allows merging multiple html files in one
interactive krona file

Note that you need to use the name args in
[`krona()`](https://adrientaudiere.github.io/MiscMetabar/reference/krona.md)
function before `merge_krona()` in order to give good name to each krona
pie in the output.

## Usage

``` r
merge_krona(files = NULL, output = "mergeKrona.html")
```

## Arguments

- files:

  (required) path to html files to merged

- output:

  path to the output file

## Value

A html file

## Details

This function is mainly a wrapper of the work of others. Please cite
[Krona](https://github.com/marbl/Krona) if you use this function.

## See also

[`krona`](https://adrientaudiere.github.io/MiscMetabar/reference/krona.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # tolower(Sys.info()[["sysname"]]) != "windows" && MiscMetabar::is_krona_installed()
if (FALSE) { # \dontrun{
data("GlobalPatterns", package = "phyloseq")
GA <- subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
krona(GA, "Number.of.sequences.html", name = "Nb_seq_GP_acidobacteria")
krona(GA, "Number.of.ASVs.html", nb_seq = FALSE, name = "Nb_asv_GP_acidobacteria")
merge_krona(c("Number.of.sequences.html", "Number.of.ASVs.html"), "mergeKrona.html")
unlink(c("Number.of.sequences.html", "Number.of.ASVs.html", "mergeKrona.html"))
} # }
}
```
