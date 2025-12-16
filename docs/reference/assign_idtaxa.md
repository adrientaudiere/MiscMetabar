# A wrapper of [`IdTaxa`](https://rdrr.io/pkg/DECIPHER/man/IdTaxa.html)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function is basically a wrapper of functions
[`DECIPHER::IdTaxa()`](https://rdrr.io/pkg/DECIPHER/man/IdTaxa.html) and
[`DECIPHER::LearnTaxa()`](https://rdrr.io/pkg/DECIPHER/man/LearnTaxa.html),
please cite the DECIPHER package if you use this function. Note that if
you want to specify parameters for the learning step you must used the
trainingSet param instead of the a fasta_for_training. The training file
can be obtain using the function
[`learn_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/learn_idtaxa.md).

It requires:

- either a physeq or seq2search object.

- either a trainingSet or a fasta_for_training

## Usage

``` r
assign_idtaxa(
  physeq,
  trainingSet = NULL,
  seq2search = NULL,
  fasta_for_training = NULL,
  behavior = "return_matrix",
  threshold = 60,
  column_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  suffix = "_idtaxa",
  nproc = 1,
  unite = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- trainingSet:

  An object of class Taxa and subclass Train compatible with the class
  of test.

- seq2search:

  A DNAStringSet object of sequences to search for. Replace the physeq
  object.

- fasta_for_training:

  A fasta file (can be gzip) to train the trainingSet using the function
  [`learn_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/learn_idtaxa.md).
  Only used if trainingSet is NULL.

  The reference database must contain taxonomic information in the
  header of each sequence in the form of a string starting with ";tax="
  and followed by a comma-separated list of up to nine taxonomic
  identifiers.

  The only exception is if `unite=TRUE`. In that case the UNITE taxonomy
  is automatically formatted.

- behavior:

  Either "return_matrix" (default), or "add_to_phyloseq":

  - "return_matrix" return a list of two objects. The first element is
    the taxonomic matrix and the second element is the raw results from
    DECIPHER::IdTaxa() function.

  - "add_to_phyloseq" return a phyloseq object with amended slot
    `@taxtable`. Only available if using physeq input and not seq2search
    input.

- threshold:

  (Int, default 60) Numeric specifying the confidence at which to
  truncate the output taxonomic classifications. Lower values of
  threshold will classify deeper into the taxonomic tree at the expense
  of accuracy, and vise-versa for higher values of threshold. See
  [`DECIPHER::IdTaxa()`](https://rdrr.io/pkg/DECIPHER/man/IdTaxa.html)
  man page.

- column_names:

  (vector of character) names for the column of the taxonomy

- suffix:

  (character) The suffix to name the new columns. Default to "\_idtaxa".

- nproc:

  (default: 1) Set to number of cpus/processors to use

- unite:

  (logical, default FALSE). If set to TRUE, the fasta_for_training file
  is formatted from UNITE format to sintax one, needed in
  fasta_for_training. Only used if trainingSet is NULL.

- verbose:

  (logical). If TRUE, print additional information.

- ...:

  Additional arguments passed on to
  [`IdTaxa`](https://rdrr.io/pkg/DECIPHER/man/IdTaxa.html)

## Value

Either a new phyloseq object with additional information in the
@tax_table slot or a list of two objects if behavior is "return_matrix"

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`DECIPHER::IdTaxa()`](https://rdrr.io/pkg/DECIPHER/man/IdTaxa.html) if
you use this function.

## See also

[`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md),
[`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md),
[`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md),
[`assign_blastn()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_blastn.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
# /!\ The value of threshold must be change for real database (recommend
#  value are between 50 and 70).

data_fungi_mini_new <- assign_idtaxa(data_fungi_mini,
  fasta_for_training = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar"
  ), threshold = 20, behavior = "add_to_phyloseq"
)

result_idtaxa <- assign_idtaxa(data_fungi_mini,
  fasta_for_training = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar"
  ), threshold = 20
)

plot(result_idtaxa$idtaxa_raw)
} # }
```
