# Add information from [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md) to the `tax_table` slot of a *phyloseq* object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Basically a wrapper of
[`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md)
with option `unique_per_seq = TRUE` and `score_filter = FALSE`.

Add the information to the taxtable

## Usage

``` r
add_blast_info(
  physeq,
  fasta_for_db,
  silent = FALSE,
  suffix = "blast_info",
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fasta_for_db:

  path to a fasta file to make the blast database

- silent:

  (logical) If true, no message are printing.

- suffix:

  (character) The suffix to name the new columns. Set the suffix to ""
  in order to remove any suffix.

- ...:

  Additional arguments passed on to
  [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md)
  function.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with more information in tax_table based on a blast on a given
database

## Author

Adrien Taudière
