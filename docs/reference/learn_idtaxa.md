# A wrapper of [`DECIPHER::LearnTaxa()`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function is basically a wrapper of functions
[`DECIPHER::LearnTaxa()`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html),
please cite the DECIPHER package if you use this function.

## Usage

``` r
learn_idtaxa(
  fasta_for_training,
  output_Rdata = NULL,
  output_path_only = FALSE,
  unite = FALSE,
  ...
)
```

## Arguments

- fasta_for_training:

  A fasta file (can be gzip) to train the trainingSet using the function
  `learn_idtaxa()`. Only used if trainingSet is NULL.

  The reference database must contain taxonomic information in the
  header of each sequence in the form of a string starting with ";tax="
  and followed by a comma-separated list of up to nine taxonomic
  identifiers.

  The only exception is if `unite=TRUE`. In that case the UNITE taxonomy
  is automatically formatted.

- output_Rdata:

  A vector naming the path to an output Rdata file. If left to NULL, no
  Rdata file is written.

- output_path_only:

  (logical, default FALSE). If TRUE, the function return only the path
  to the output_Rdata file. Note that output_Rdata must be set.

- unite:

  (logical, default FALSE). If set to TRUE, the fasta_for_training file
  is formatted from UNITE format to sintax one, needed in
  fasta_for_training. Only used if trainingSet is NULL.

- ...:

  Additional arguments passed on to
  [`DECIPHER::LearnTaxa()`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html)

## Value

Either a Taxa Train object (see
[`DECIPHER::LearnTaxa()`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html))
or, if output_path_only is TRUE, a vector indicating the path to the
output training object.

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`DECIPHER::LearnTaxa()`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html)
if you use this function.

## See also

[`assign_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_idtaxa.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
training_mini_UNITE_fungi <-
  learn_idtaxa(fasta_for_training = system.file("extdata",
    "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar"
  ))
plot(training_mini_UNITE_fungi)

training_100sp_UNITE <-
  learn_idtaxa(
    fasta_for_training = system.file("extdata",
      "100_sp_UNITE_sh_general_release_dynamic.fasta",
      package = "MiscMetabar"
    ),
    unite = TRUE
  )

plot(training_100sp_UNITE)
} # }
```
