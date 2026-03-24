# Format a fasta database in dada2 format

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

First format in sintax format and then in dada2 format

## Usage

``` r
format2dada2(
  fasta_db = NULL,
  taxnames = NULL,
  output_path = NULL,
  from_sintax = TRUE,
  pattern_to_remove = NULL,
  ...
)
```

## Arguments

- fasta_db:

  A link to a fasta files

- taxnames:

  A list of names to format. You must specify either fasta_db OR
  taxnames, not both.

- output_path:

  (optional) A path to an output fasta files. Only used if fasta_db is
  set.

- from_sintax:

  (logical, default FALSE) Is the original fasta file in sintax format?

- pattern_to_remove:

  (a regular expression) Define a pattern to remove. For example,
  pattern_to_remove = "\\rep.\*" remove all character after '\|rep' to
  force
  [`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html)
  to not use the database as a Unite-formated database

- ...:

  Additional arguments passed on to
  [`format2sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2sintax.md)
  function

## Value

Either an object of class DNAStringSet or a vector of reformated names

## See also

[`format2dada2_species()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2dada2_species.md),
[`format2sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2sintax.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
f <- system.file("extdata", "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar"
)
format2dada2(fasta_db = f, from_sintax = FALSE)
#> DNAStringSet object of length 5000:
#>        width seq                                            names               
#>    [1]   492 GGAAGGATCATTAACACTAATC...GCTTAGTCACACAAACTTTAC k__Fungi;p__Fungi...
#>    [2]   451 CGAGTTTACCACTCCAACCCCT...GTTAAACCCCAACTTCTGAAT k__Fungi;p__Ascom...
#>    [3]   457 CCGAGTTTACAACTCCCAAGCC...CAGTTAACCCTAAGTTTTGAT k__Fungi;p__Fungi...
#>    [4]   490 AGAAGGGGGAAGAGCAATTCTA...GTTCAACCCCAACTTCTGAAT k__Fungi;p__Ascom...
#>    [5]   464 CCAGTTGCAACATACCACTTGT...GTTAAACCCCAACTTCTGAAT k__Fungi;p__Ascom...
#>    ...   ... ...
#> [4996]   622 CATTATTGAACAAAAAAATGGT...GATTTGGAGAACCAAGATAAA k__Fungi;p__Basid...
#> [4997]   609 ATGAATTGTAACGGGGGTTGAT...CTGATCTGTTCCACTTTCAAC k__Fungi;p__Basid...
#> [4998]   615 ATGGAATTGTAACAGGAGCTGA...AATTCTGTTCAACTTCCACCC k__Fungi;p__Fungi...
#> [4999]   606 ATGAATTGTAACAGGGGCTGAT...TGATCTGTTCCATTTTTCAAC k__Fungi;p__Basid...
#> [5000]   519 TTGGAGTCTGGTTGTTGCTGGC...ATTTTCATTTCATCCTTGCCA k__Fungi;p__Basid...
# }
```
