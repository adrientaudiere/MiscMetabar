# Subsample a fastq file copying the n_seq first sequences in a given folder

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful to test a pipeline on small fastq files.

## Usage

``` r
subsample_fastq(fastq_files, folder_output = "subsample", nb_seq = 1000)
```

## Arguments

- fastq_files:

  The path to one fastq file or a list of fastq files (see examples)

- folder_output:

  The path to a folder for output files

- nb_seq:

  (int; default 1000) : Number of sequences kept (every sequence spread
  across 4 lines)

## Value

Nothing, create subsampled fastq files in a folder

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
ex_file <- system.file("extdata", "ex_R1_001.fastq.gz",
  package = "MiscMetabar",
  mustWork = TRUE
)
subsample_fastq(ex_file, paste0(tempdir(), "/output_fastq"))
subsample_fastq(list_fastq_files(system.file("extdata", package = "MiscMetabar")),
  paste0(tempdir(), "/output_fastq"),
  n = 10
)
unlink(paste0(tempdir(), "/output_fastq"), recursive = TRUE)
# }
```
