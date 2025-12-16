# Detect for chimera taxa using [vsearch](https://github.com/torognes/vsearch)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the VSEARCH software.

## Usage

``` r
chimera_detection_vs(
  seq2search,
  nb_seq,
  vsearchpath = "vsearch",
  abskew = 2,
  min_seq_length = 100,
  vsearch_args = "--fasta_width 0",
  keep_temporary_files = FALSE
)
```

## Arguments

- seq2search:

  (required) a list of DNA sequences coercible by function
  [`Biostrings::DNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)

- nb_seq:

  (required) a numeric vector giving the number of sequences for each
  DNA sequences

- vsearchpath:

  (default: "vsearch") path to vsearch

- abskew:

  (int, default 2) The abundance skew is used to distinguish in a three
  way alignment which sequence is the chimera and which are the parents.
  The assumption is that chimeras appear later in the PCR amplification
  process and are therefore less abundant than their parents. The
  default value is 2.0, which means that the parents should be at least
  2 times more abundant than their chimera. Any positive value equal or
  greater than 1.0 can be used.

- min_seq_length:

  (int, default 100)) Minimum length of sequences to be part of the
  analysis

- vsearch_args:

  (default "–fasta_width 0") A list of other args for vsearch command

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files ?

  - non_chimeras.fasta

  - chimeras.fasta

  - borderline.fasta

## Value

A list of 3 including non-chimera taxa (`$non_chimera`), chimera taxa
(`$chimera`) and bordeline taxa (`$borderline`)

## Details

This function is mainly a wrapper of the work of others. Please make
[vsearch](https://github.com/torognes/vsearch).

## See also

[`chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/reference/chimera_removal_vs.md),
[`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
chimera_detection_vs(
  seq2search = data_fungi@refseq,
  nb_seq = taxa_sums(data_fungi)
)
#> Filtering for sequences under 100 bp remove a total of 0 ( 0 %) unique sequences for a total of 0 sequences removed ( 0 %)
#> $non_chimera
#> AAStringSet object of length 1054:
#>        width seq                                            names               
#>    [1]   312 AAATGCGATAAGTAATGTGAAT...TAGGAATACCCGCTGAACTTA Taxa1;size=92884
#>    [2]   301 AAATGCGATAAGTAATGTGAAT...TAGGAATACCCGCTGAACTTA Taxa2;size=53538
#>    [3]   349 AAATGCGATAAGTAATGTGAAT...TGGGACTACCCGCTGAACTTA Taxa3;size=47410
#>    [4]   357 AAATGCGATAAGTAATGTGAAT...TGGGACTACCCGCTGAACTTA Taxa4;size=46857
#>    [5]   300 AAATGCGATAAGTAATGTGAAT...TAGGAATACCCGCTGAACTTA Taxa5;size=41082
#>    ...   ... ...
#> [1050]   260 AAACGCGAAAAGTGTTATGATG...AAGATCACCCGCTGAACTTAA Taxa1420;size=2
#> [1051]   365 AAATGCGATAAGTAATGTGAAT...TAGGACTACCCGCTGAACTTA Taxa602;size=2
#> [1052]   344 AAATGCGATAAGTAATGTGAAT...TAGGAATACGCGCTGAACTTA Taxa1142;size=1
#> [1053]   290 GAAATGCGATAAGTAATGTGAA...TAGGGATACCCGCTGAACTTA Taxa246;size=1
#> [1054]   318 GAAATGCGATACGTAATGTGAA...AGGGATACCCGCTGAACTTAA Taxa412;size=1
#> 
#> $chimera
#> AAStringSet object of length 240:
#>       width seq                                             names               
#>   [1]   341 GAAATGCGATAAGTAATGTGAA...GTAGGATTACCCGCTGAACTTA Taxa136;size=2743
#>   [2]   307 AAATGCGATAAGTAATGTGAAT...GTAGGGATACCCGCTGAACTTA Taxa163;size=2028
#>   [3]   339 AAATGCGATAAGTAATGTGAAT...GTAGGATTACCCGCTGAACTTA Taxa206;size=1471
#>   [4]   312 AAATGCGATAAGTAATGTGAAT...GTAGGAATACCCGCTGAACTTA Taxa286;size=875
#>   [5]   306 AAATGCGAAAAGTAGTGTGAAT...GTAGGGATACCCGCTGAACTTA Taxa294;size=832
#>   ...   ... ...
#> [236]   303 GAAATGCGATACTTGGTGTGAA...GTGGGACTACCCGCTGAACTTA Taxa1390;size=29
#> [237]   293 GAACTACGATAAGTAATGTGAA...TAGGGATACCCGCTGAACTTAA Taxa1395;size=25
#> [238]   301 AAATGCGAAAAGTAGTGTGAAT...TAGGGATACCCGCTGAACTTAA Taxa1407;size=18
#> [239]   382 AAATGCGATAAATAATATGAAT...GCAAGATTACCCGCTGAACTTA Taxa1411;size=12
#> [240]   287 GAAATGCGATAAGTAATGTGAA...TAGGGCTACCCGCTGAACTTAA Taxa789;size=2
#> 
#> $borderline
#> AAStringSet object of length 126:
#>       width seq                                             names               
#>   [1]   298 AAATGCGATAAGTAATGTGAAT...GTAGGGATACCCGCTGAACTTA Taxa123;size=3164
#>   [2]   300 AAATGCGATAAGTAATGTGAAT...GTAGGGATACCCGCTGAACTTA Taxa139;size=2712
#>   [3]   304 AAATGCGATAAGTAATGTGAAT...GTAGGACTACCCGCTGAACTTA Taxa141;size=2611
#>   [4]   295 AAATGCGATAAGTAATGTGAAT...GTAGGGATACCCGCTGAACTTA Taxa153;size=2294
#>   [5]   302 AAATGCGATAAGTAATGTGAAT...GTAGGGATACCCGCTGAACTTA Taxa174;size=1851
#>   ...   ... ...
#> [122]   330 AAACGCGATAGGTAATGTGAAT...GTAGGACTACCCGCTGAACTTA Taxa1292;size=62
#> [123]   332 AAATGCGATAAGTAATGTGAAT...GTAGGAATACCCGCTGAACTTA Taxa1353;size=48
#> [124]   300 GAAATGCGATAAGTAATGCGAA...GTAGGGATACCCGCTGAACTTA Taxa1384;size=31
#> [125]   302 GAAATGCGATACGTAATGTGAA...TAGGAATACCCGCTGAACTTAA Taxa1388;size=29
#> [126]   316 GAAATGCGATAAATAATATGAA...ACAAGATTACCCGCTGAACTTA Taxa1399;size=22
#> 
# }
```
