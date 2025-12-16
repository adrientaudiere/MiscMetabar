# Filter undesirable taxa using blast against a custom database.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the blast software.

## Usage

``` r
filter_asv_blast(
  physeq,
  fasta_for_db = NULL,
  database = NULL,
  clean_pq = TRUE,
  add_info_to_taxtable = TRUE,
  id_filter = 90,
  bit_score_filter = 50,
  min_cover_filter = 50,
  e_value_filter = 1e-30,
  ...
)

filter_taxa_blast(
  physeq,
  fasta_for_db = NULL,
  database = NULL,
  clean_pq = TRUE,
  add_info_to_taxtable = TRUE,
  id_filter = 90,
  bit_score_filter = 50,
  min_cover_filter = 50,
  e_value_filter = 1e-30,
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

- database:

  path to a blast database

- clean_pq:

  (logical) If set to TRUE, empty samples and empty taxa (ASV, OTU) are
  discarded after filtering.

- add_info_to_taxtable:

  (logical, default TRUE) Does the blast information are added to the
  taxtable ?

- id_filter:

  (default: 90) cut of in identity percent to keep result

- bit_score_filter:

  (default: 50) cut of in bit score to keep result The higher the
  bit-score, the better the sequence similarity. The bit-score is the
  requires size of a sequence database in which the current match could
  be found just by chance. The bit-score is a log2 scaled and normalized
  raw-score. Each increase by one doubles the required database size
  (2bit-score).

- min_cover_filter:

  (default: 50) cut of in query cover (%) to keep result

- e_value_filter:

  (default: 1e-30) cut of in e-value (%) to keep result The BLAST
  E-value is the number of expected hits of similar quality (score) that
  could be found just by chance.

- ...:

  Additional arguments passed on
  to[`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md)
  function. See
  [`?blast_pq`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md).
  Note that params `unique_per_seq` must be left to TRUE and
  `score_filter` must be left to FALSE.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object.
