# Dereplicate sequences in a fastx file

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wrapper of
[`Rsearch::vs_fastx_uniques()`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_uniques.html)
(vsearch `--fastx_uniques`). Dereplicate the reads of a fasta or fastq
file: identical sequences are merged and their abundances recorded as
`;size=N` annotations in the output headers.

## Usage

``` r
vs_fastx_uniques(
  fastx_input,
  output_format = "fasta",
  fastx_output = NULL,
  minuniquesize = 1,
  strand = "plus",
  sizein = TRUE,
  sizeout = TRUE,
  relabel = NULL,
  relabel_sha1 = FALSE,
  fastq_qout_max = FALSE,
  fasta_width = 0,
  sample = NULL,
  vsearch_options = NULL,
  tmpdir = NULL
)
```

## Arguments

- fastx_input:

  (required) Path to a fasta or fastq file (optionally gzipped).

- output_format:

  (default `"fasta"`) Output format: `"fasta"` or `"fastq"`.

- fastx_output:

  Path to the output file. A temporary file is used when `NULL`.

- minuniquesize:

  (integer, default 1) Discard unique sequences with an abundance lower
  than this value.

- strand:

  (default `"plus"`) Compare sequences on the `"plus"` strand only or on
  `"both"` strands.

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

- fastq_qout_max:

  (logical, default FALSE) Report the maximum quality score in fastq
  output.

- fasta_width:

  (integer, default 0) Width of fasta sequence lines (0 for no
  wrapping).

- sample:

  Optional sample name annotation added to output headers.

- vsearch_options:

  Additional vsearch command-line options.

- tmpdir:

  Directory for temporary files.

## Value

`NULL` (invisibly). The dereplicated sequences are written to
`fastx_output`.

## See also

[`Rsearch::vs_fastx_uniques()`](https://cassandrahjo.github.io/Rsearch/reference/vs_fastx_uniques.html),
[`vs_uchime_ref()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/vs_uchime_ref.md),
[`cluster_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cluster_reads.md)

## Author

Adrien Taudière

## Examples

``` r
vs_fastx_uniques(
  system.file("extdata", "ex_little.fasta", package = "MiscMetabar"),
  fastx_output = tempfile(fileext = ".fasta")
)
```
