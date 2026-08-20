# Reference-based chimera detection with vsearch

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wrapper of
[`Rsearch::vs_uchime_ref()`](https://cassandrahjo.github.io/Rsearch/reference/vs_uchime_ref.html)
(vsearch `--uchime_ref`). Detect chimeric sequences in a fasta file by
comparison against a reference database. Complements
[`chimera_detection_vs()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/chimera_detection_vs.md),
which performs *de novo* detection without a reference database.

## Usage

``` r
vs_uchime_ref(
  fasta_input,
  database,
  nonchimeras = NULL,
  chimeras = NULL,
  sizein = TRUE,
  sizeout = TRUE,
  relabel = NULL,
  relabel_sha1 = FALSE,
  fasta_width = 0,
  sample = NULL,
  log_file = NULL,
  threads = 1,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fasta_input:

  (required) Path to the query fasta file.

- database:

  (required) Path to the reference database fasta file.

- nonchimeras:

  Path to the output fasta file of non-chimeric sequences. A temporary
  file is used when `NULL`.

- chimeras:

  Path to the output fasta file of chimeric sequences. A temporary file
  is used when `NULL`.

- sizein:

  (logical, default TRUE) Take `;size=N` annotations of the input into
  account.

- sizeout:

  (logical, default TRUE) Add `;size=N` annotations to the output.

- relabel:

  Relabel output sequences with this prefix followed by a sequential
  number.

- relabel_sha1:

  (logical, default FALSE) Relabel output sequences with their SHA1
  digest.

- fasta_width:

  (integer, default 0) Width of fasta sequence lines (0 for no
  wrapping).

- sample:

  Optional sample name annotation added to output headers.

- log_file:

  Path to a log file.

- threads:

  (integer, default 1) Number of threads.

- vsearch_options:

  Additional vsearch command-line options.

- tmpdir:

  Directory for temporary files.

## Value

`NULL` (invisibly). The non-chimeric and chimeric sequences are written
to the `nonchimeras` and `chimeras` files, respectively.

## See also

[`Rsearch::vs_uchime_ref()`](https://cassandrahjo.github.io/Rsearch/reference/vs_uchime_ref.html),
[`chimera_detection_vs()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/chimera_detection_vs.md),
[`chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/chimera_removal_vs.md)

## Author

Adrien Taudière

## Examples

``` r
vs_uchime_ref(
  system.file("extdata", "ex.fasta", package = "MiscMetabar"),
  database = system.file(
    "extdata",
    "100_sp_UNITE_sh_general_release_dynamic.fasta",
    package = "MiscMetabar"
  )
)
#> # A tibble: 3 × 2
#>   Header                 Sequence                                               
#>   <chr>                  <chr>                                                  
#> 1 Eutypa_spinosa;size=1  TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTAC…
#> 2 Eutypa_spinosa2;size=1 TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTAC…
#> 3 Eutypa_spinosa3;size=1 TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCT…
```
