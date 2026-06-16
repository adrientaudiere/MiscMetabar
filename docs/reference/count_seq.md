# Count sequences in fasta or fastq file

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

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
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex.fasta
#> [1] 3
count_seq(
  folder_path = system.file("extdata", package = "MiscMetabar"),
  pattern = "*.fasta"
)
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/100_sp_UNITE_sh_general_release_dynamic.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/100_sp_UNITE_sh_general_release_dynamic.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/100_sp_UNITE_sh_general_release_dynamic_dadaSpecies.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/100_sp_UNITE_sh_general_release_dynamic_dadaSpecies.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/100_sp_UNITE_sh_general_release_dynamic_sintax.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/100_sp_UNITE_sh_general_release_dynamic_sintax.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex_R1_001.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex_R1_001.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex_R1_002.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex_R1_002.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex_little.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/ex_little.fasta
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/mini_UNITE_fungi.fasta.gz
#> Warning: There is more than one '.' inside your file path: /home/adrien/R/x86_64-pc-linux-gnu-library/4.6/MiscMetabar/extdata/mini_UNITE_fungi.fasta.gz
#> [1]  100  100  100    3   63   51    2 5000
```
