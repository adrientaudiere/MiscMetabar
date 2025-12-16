# Retrieve the FUNGuild database

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The original function and documentation was written by Brendan Furneaux
in the [FUNGuildR](https://github.com/brendanf/FUNGuildR/) package.

Please cite this publication
([doi:10.1016/j.funeco.2015.06.006](https://doi.org/10.1016/j.funeco.2015.06.006)
).

## Usage

``` r
get_funguild_db(db_url = "http://www.stbates.org/funguild_db_2.php")
```

## Arguments

- db_url:

  a length 1 character string giving the URL to retrieve the database
  from

## Value

a [`tibble::tibble`](https://tibble.tidyverse.org/reference/tibble.html)
containing the database, which can be passed to the `db` argument of
[`funguild_assign()`](https://adrientaudiere.github.io/MiscMetabar/reference/funguild_assign.md)

## References

Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J, Schilling
JS, Kennedy PG. 2016. *FUNGuild: An open annotation tool for parsing
fungal community datasets by ecological guild*. Fungal Ecology
20:241-248.

## Author

Brendan Furneaux (orcid:
[0000-0003-3522-7363](https://orcid.org/0000-0003-3522-7363)), modified
by Adrien Taudière
