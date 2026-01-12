# Count sequences in fasta or fastq file

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use grep to count the number of line with only one '+' (fastq, fastq.gz)
or lines starting with a '\>' (fasta) to count sequences.

## Usage

``` r
count_seq(file_path = NULL, folder_path = NULL, pattern = NULL)
```

## Arguments

- file_path:

  The path to a fasta, fastq or fastq.gz file

- folder_path:

  The path to a folder with fasta, fastq or fastq.gz files

- pattern:

  A pattern to filter files in a folder. E.g. *R2*

## Value

the number of sequences

## Author

Adrien Taudière

## Examples

``` r
count_seq(file_path = system.file(
  "extdata",
  "ex.fasta",
  package = "MiscMetabar",
  mustWork = TRUE
))
#> [1] 3
count_seq(
  folder_path = system.file("extdata", package = "MiscMetabar"),
  pattern = "*.fasta"
)
#> Warning: There is more than one '.' inside your file path: /tmp/RtmpV46sfz/temp_libpath1e51304f43f9/MiscMetabar/extdata/mini_UNITE_fungi.fasta.gz
#> Warning: There is more than one '.' inside your file path: /tmp/RtmpV46sfz/temp_libpath1e51304f43f9/MiscMetabar/extdata/mini_UNITE_fungi.fasta.gz
#> [1]  100  100  100    3   63   51    2 5000
```
