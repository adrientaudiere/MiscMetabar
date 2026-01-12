# List fastq files

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Useful for targets bioinformatic pipeline.

## Usage

``` r
list_fastq_files(
  path,
  paired_end = TRUE,
  pattern = "fastq",
  pattern_R1 = "_R1_",
  pattern_R2 = "_R2_",
  nb_files = Inf
)
```

## Arguments

- path:

  path to files (required)

- paired_end:

  do you have paired_end files? (default TRUE)

- pattern:

  a pattern to filter files (passed on to list.files function).

- pattern_R1:

  a pattern to filter R1 files (default "*R1*")

- pattern_R2:

  a pattern to filter R2 files (default "*R2*")

- nb_files:

  the number of fastq files to list (default FALSE)

## Value

a list of one (single end) or two (paired end) list of files files are
sorted by names (default behavior of
[`list.files()`](https://rdrr.io/r/base/list.files.html))

## Author

Adrien Taudière

## Examples

``` r
list_fastq_files(system.file("extdata", package = "MiscMetabar"))
#> $fnfs
#> [1] "/tmp/RtmpV46sfz/temp_libpath1e51304f43f9/MiscMetabar/extdata/ex_R1_001.fastq.gz"
#> 
#> $fnrs
#> [1] "/tmp/RtmpV46sfz/temp_libpath1e51304f43f9/MiscMetabar/extdata/ex_R2_001.fastq.gz"
#> 
list_fastq_files(system.file("extdata", package = "MiscMetabar"),
  paired_end = FALSE, pattern_R1 = ""
)
#> $fnfs
#> [1] "/tmp/RtmpV46sfz/temp_libpath1e51304f43f9/MiscMetabar/extdata/ex.fastq"          
#> [2] "/tmp/RtmpV46sfz/temp_libpath1e51304f43f9/MiscMetabar/extdata/ex_R1_001.fastq.gz"
#> [3] "/tmp/RtmpV46sfz/temp_libpath1e51304f43f9/MiscMetabar/extdata/ex_R2_001.fastq.gz"
#> 
```
