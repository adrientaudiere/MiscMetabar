# Add information about Guild for FUNGI the FUNGuild databse

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Please cite Nguyen et al. 2016
([doi:10.1016/j.funeco.2015.06.006](https://doi.org/10.1016/j.funeco.2015.06.006)
)

## Usage

``` r
add_funguild_info(
  physeq,
  taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  db_url = "http://www.stbates.org/funguild_db_2.php"
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- taxLevels:

  Name of the 7 columns in tax_table required by funguild

- db_url:

  a length 1 character string giving the URL to retrieve the database
  from

## Value

A new object of class `physeq` with Guild information added to
`tax_table` slot

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to `FUNGuildR` package and the associate publication
([doi:10.1016/j.funeco.2015.06.006](https://doi.org/10.1016/j.funeco.2015.06.006)
) if you use this function.

## See also

[`plot_guild_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_guild_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# to avoid bug in CRAN when internet is not available
if (requireNamespace("httr")) {
  d_fung_mini <- add_funguild_info(data_fungi_mini,
    taxLevels = c(
      "Domain",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )
  )
  sort(table(d_fung_mini@tax_table[, "guild"]), decreasing = TRUE)
}
} # }
```
