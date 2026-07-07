# Convert FASTQ file(s) to FASTA file(s)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Convert one or several FASTQ files (optionally gzip-compressed) to FASTA
format. For each input file, the `@` sign starting a FASTQ sequence
header is replaced by the `>` sign required by the FASTA format, and the
quality lines are dropped. Output files keep the input base name with
the extension changed to `.fasta`.

Inspired by <https://github.com/ambuvjyn/baseq>.

## Usage

``` r
fastq_to_fasta(fastq_files, output_folder = tempdir(), force = FALSE)
```

## Arguments

- fastq_files:

  (character, required) One or more paths to FASTQ files. The `.fastq`,
  `.fq`, `.fastq.gz` and `.fq.gz` extensions are recognized
  (case-insensitive).

- output_folder:

  (character, default:
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html)) Folder where the
  FASTA files are written. Created if it does not exist.

- force:

  (logical, default: FALSE) If FALSE, existing output files are not
  overwritten and an error is thrown instead.

## Value

A character vector of the paths to the output FASTA files (invisibly
named by the input files).

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
fq <- system.file("extdata", "ex.fastq", package = "MiscMetabar")
out <- fastq_to_fasta(fq, output_folder = tempdir())
#> ✔ Wrote 1 FASTA file to /tmp/Rtmp2sWBuj.
readLines(out, n = 2)
#> [1] ">M03493:483:000000000-KW735:1:1101:10710:1281 1:N:0:17"                                                                                                                                                                                                                                                       
#> [2] "CTTGTTGGAGGGCAAGTCTGGTGCCCGCAGCCGCGGTAATTCCAGCTCCAATAGCGNNTNNTAAAGTTGTTGCAGTTAAAAAGCTCGTAGTTGAATTTCGGGATNNATGGGTTGGTCACGCCCTCGGTATGTACTGATTGCATTGGTGCTCACCTTCTAAGGAACCATGATGTCATTCGTTTGATGTGGTTGGGAATGAGGACTGTTATCTTGAAAAAATTGGAGGGTTTAAAGCAGTCTAACGCTTTAATACATTATCATGGAATAATGAAATCGGACAGCCGATTATATTTTGTTGGCG"
# }
```
