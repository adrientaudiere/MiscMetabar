# Assign Guilds to Organisms Based on Taxonomic Classification

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The original function and documentation was written by Brendan Furneaux
in the [FUNGuildR](https://github.com/brendanf/FUNGuildR/) package.

These functions have identical behavior if supplied with a database;
however they download the database corresponding to their name by
default.

Taxa present in the database are matched to the taxa present in the
supplied `otu_table` by exact name. In the case of multiple matches, the
lowest (most specific) rank is chosen. No attempt is made to check or
correct the classification in `otu_table$Taxonomy`.

## Usage

``` r
funguild_assign(
  otu_table,
  db_url = NULL,
  db_funguild = NULL,
  tax_col = "Taxonomy"
)
```

## Arguments

- otu_table:

  A `data.frame` with a `character` column named "`Taxonomy`" (or
  another name as specified in `tax_col`), as well as any other columns.
  Each entry in "`otu_table$Taxonomy`" should be a comma-, colon-,
  underscore-, or semicolon-delimited classification of an organism.
  Rank indicators as given by Sintax ("`k:`", "`p:`"...) or Unite
  ("`k__`, "`p__`", ...) are also allowed. A `character` vector,
  representing only the taxonomic classification, is also accepted.

- db_url:

  a length 1 character string giving the URL to retrieve the database
  from

- db_funguild:

  A `data.frame` representing the FUNGuild as returned by
  [`get_funguild_db()`](https://adrientaudiere.github.io/MiscMetabar/reference/get_funguild_db.md)
  If not supplied, the default database will be downloaded.

- tax_col:

  A `character` string, optionally giving an alternate column name in
  `otu_table` to use instead of `otu_table$Taxonomy`.

## Value

A [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing all columns of `otu_table`, plus relevant columns of
information from the FUNGuild

## References

Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J, Schilling
JS, Kennedy PG. 2016. *FUNGuild: An open annotation tool for parsing
fungal community datasets by ecological guild*. Fungal Ecology
20:241-248.

## Author

Brendan Furneaux (orcid:
[0000-0003-3522-7363](https://orcid.org/0000-0003-3522-7363)), modified
by Adrien Taudière
