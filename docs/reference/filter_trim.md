# A wrapper of the function [`dada2::filterAndTrim()`](https://rdrr.io/pkg/dada2/man/filterAndTrim.html) to use in [targets](https://books.ropensci.org/targets/) pipeline

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function filter and trim (with parameters passed on to
[`dada2::filterAndTrim()`](https://rdrr.io/pkg/dada2/man/filterAndTrim.html)
function) forward sequences or paired end sequence if 'rev' parameter is
set. It return the list of files to subsequent analysis in a targets
pipeline.

## Usage

``` r
filter_trim(
  fw = NULL,
  rev = NULL,
  output_fw = file.path(paste(getwd(), "/output/filterAndTrim_fwd", sep = "")),
  output_rev = file.path(paste(getwd(), "/output/filterAndTrim_rev", sep = "")),
  return_a_vector = FALSE,
  ...
)
```

## Arguments

- fw:

  (required) a list of forward fastq files

- rev:

  a list of reverse fastq files for paired end trimming

- output_fw:

  Path to output folder for forward files. By default, this function
  will create a folder "output/filterAndTrim_fwd" in the current working
  directory.

- output_rev:

  Path to output folder for reverse files. By default, this function
  will create a folder "output/filterAndTrim_fwd" in the current working
  directory.

- return_a_vector:

  (logical, default FALSE) If true, the return is a vector of path
  (usefull when used with targets::tar_targets(..., format="file"))

- ...:

  Other parameters passed on to
  [`dada2::filterAndTrim()`](https://rdrr.io/pkg/dada2/man/filterAndTrim.html)
  function.

## Value

A list of files. If rev is set, will return a list of two lists. The
first list is a list of forward files, and the second one is a list of
reverse files.

## See also

[`dada2::filterAndTrim()`](https://rdrr.io/pkg/dada2/man/filterAndTrim.html)

## Author

Adrien Taudière

## Examples

``` r
testFastqs_fw <- c(
  system.file("extdata", "sam1F.fastq.gz", package = "dada2"),
  system.file("extdata", "sam2F.fastq.gz", package = "dada2")
)
testFastqs_rev <- c(
  system.file("extdata", "sam1R.fastq.gz", package = "dada2"),
  system.file("extdata", "sam2R.fastq.gz", package = "dada2")
)

filt_fastq_fw <- filter_trim(testFastqs_fw, output_fw = tempdir())
derep_fw <- derepFastq(filt_fastq_fw[1])
derep_fw
#> $sam1F.fastq.gz
#> derep-class: R object describing dereplicated sequencing reads
#> $uniques: 1500 reads in 896 unique sequences
#>   Sequence lengths: min=250, median=250, max=250
#> $quals: Quality matrix dimension:  896 250
#>   Consensus quality scores: min=7, median=36, max=38
#> $map: Map from reads to unique sequences:  4 155 627 265 5 ...
#> 
#> $sam2F.fastq.gz
#> derep-class: R object describing dereplicated sequencing reads
#> $uniques: 1500 reads in 909 unique sequences
#>   Sequence lengths: min=250, median=250, max=250
#> $quals: Quality matrix dimension:  909 250
#>   Consensus quality scores: min=7, median=36, max=38
#> $map: Map from reads to unique sequences:  890 14 2 246 797 ...
#> 

if (FALSE) { # \dontrun{
filt_fastq_pe <- filter_trim(testFastqs_fw,
  testFastqs_rev,
  output_fw = paste0(tempdir(), "/", "fw"),
  output_rev = paste0(tempdir(), "rev")
)
derep_fw_pe <- derepFastq(filt_fastq_pe[[1]])
derep_rv_pe <- derepFastq(filt_fastq_pe[[2]])
derep_fw_pe
derep_rv_pe
} # }
```
