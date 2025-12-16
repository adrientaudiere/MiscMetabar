# Assign taxonomy with dada2 using 2 steps assignTaxonomy and assignSpecies

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mainly a wrapper of
[`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html)
and
[`dada2::assignSpecies()`](https://rdrr.io/pkg/dada2/man/assignSpecies.html)

## Usage

``` r
assign_dada2(
  physeq = NULL,
  ref_fasta = NULL,
  seq2search = NULL,
  min_bootstrap = 0.5,
  tryRC = FALSE,
  taxa_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    "taxId"),
  use_assignSpecies = TRUE,
  trunc_absent_ranks = FALSE,
  nproc = 1,
  suffix = "",
  verbose = TRUE,
  seq_at_one_time = 2000,
  allowMultiple = FALSE,
  from_sintax = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ref_fasta:

  (required) A link to a database in fasta.

- seq2search:

  A DNAStringSet object of sequences to search for. Replace the physeq
  object.

- min_bootstrap:

  (Float \[0:1\], default 0.5), See
  [`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html)

- tryRC:

  See
  [`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html)

- taxa_ranks:

  (vector of character) names for the column of the taxonomy

- use_assignSpecies:

  (logical, default TRUE) Do the Species rank is obtained using
  [`dada2::assignSpecies()`](https://rdrr.io/pkg/dada2/man/assignSpecies.html)
  ?

- trunc_absent_ranks:

  (logical, default FALSE) Do ranks present in taxa_ranks but not
  present in the database are removed ?

- nproc:

  (Float \[0:1\], default 0.5)

- suffix:

  (character) The suffix to name the new columns. Default to "\_idtaxa".

- verbose:

  (logical). If TRUE, print additional information.

- seq_at_one_time:

  How many sequences are treated at one time. See param `n` in
  [`dada2::assignSpecies()`](https://rdrr.io/pkg/dada2/man/assignSpecies.html)

- allowMultiple:

  (logical, default FALSE). Unchanged from
  [`dada2::assignSpecies()`](https://rdrr.io/pkg/dada2/man/assignSpecies.html).
  Defines the behavior when multiple exact matches against different
  species are returned. By default only unambiguous identifications are
  return. If TRUE, a concatenated string of all exactly matched species
  is returned. If an integer is provided, multiple identifications up to
  that many are returned as a concatenated string.

- from_sintax:

  (logical, default FALSE). Set to TRUE if the ref_fasta database is in
  sintax format. See
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md)
  for more information about the sintax format.

## Value

Either a an object of class phyloseq (if `physeq` is not NULL), or a
taxonomic table if `seq2search` is used in place of `physeq`

## Examples

``` r
if (FALSE) { # \dontrun{
data_fungi_mini2 <- assign_dada2(data_fungi_mini,
  ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar"
  ), suffix = "_dada2",
  from_sintax = TRUE
)
} # }
```
