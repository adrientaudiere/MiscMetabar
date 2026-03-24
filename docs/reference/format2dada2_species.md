# Format a fasta database in dada2 format for Species assignment

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

First format in sintax format and then in dada2 format

## Usage

``` r
format2dada2_species(
  fasta_db = NULL,
  taxnames = NULL,
  from_sintax = FALSE,
  output_path = NULL,
  ...
)
```

## Arguments

- fasta_db:

  A link to a fasta files

- taxnames:

  A list of names to format. You must specify either fasta_db OR
  taxnames, not both.

- from_sintax:

  (logical, default FALSE) Is the original fasta file in sintax format?

- output_path:

  (optional) A path to an output fasta files. Only used if fasta_db is
  set.

- ...:

  Additional arguments passed on to
  [`format2sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2sintax.md)
  function

## Value

Either an object of class DNAStringSet or a vector of reformated names

## See also

`format2dada2_species()`,
[`format2sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/format2sintax.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
f <- system.file("extdata", "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar"
)
format2dada2_species(fasta_db = f)
#> DNAStringSet object of length 5000:
#>        width seq                                            names               
#>    [1]   492 GGAAGGATCATTAACACTAATC...GCTTAGTCACACAAACTTTAC Fungi_sp_LR994329...
#>    [2]   451 CGAGTTTACCACTCCAACCCCT...GTTAAACCCCAACTTCTGAAT Fusarium_oxysporu...
#>    [3]   457 CCGAGTTTACAACTCCCAAGCC...CAGTTAACCCTAAGTTTTGAT Fungi_sp_KP861831...
#>    [4]   490 AGAAGGGGGAAGAGCAATTCTA...GTTCAACCCCAACTTCTGAAT Fusarium_prolifer...
#>    [5]   464 CCAGTTGCAACATACCACTTGT...GTTAAACCCCAACTTCTGAAT Fusarium_oxysporu...
#>    ...   ... ...
#> [4996]   622 CATTATTGAACAAAAAAATGGT...GATTTGGAGAACCAAGATAAA Amanita_pseudoscu...
#> [4997]   609 ATGAATTGTAACGGGGGTTGAT...CTGATCTGTTCCACTTTCAAC Clavulina_sp_UDB0...
#> [4998]   615 ATGGAATTGTAACAGGAGCTGA...AATTCTGTTCAACTTCCACCC Fungi_sp_LC482542...
#> [4999]   606 ATGAATTGTAACAGGGGCTGAT...TGATCTGTTCCATTTTTCAAC Clavulina_sp_UDB0...
#> [5000]   519 TTGGAGTCTGGTTGTTGCTGGC...ATTTTCATTTCATCCTTGCCA Termitomyces_sp_M...
# }
```
