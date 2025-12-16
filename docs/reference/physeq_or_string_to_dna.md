# Return a DNAStringSet object from either a character vector of DNA sequences or the `refseq` slot of a phyloseq-class object

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Internally used in
[`vsearch_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/vsearch_clustering.md),
[`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/swarm_clustering.md)
and
[`postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md).

## Usage

``` r
physeq_or_string_to_dna(physeq = NULL, dna_seq = NULL)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- dna_seq:

  You may directly use a character vector of DNA sequences in place of
  physeq args. When physeq is set, dna sequences take the value of
  `physeq@refseq`

## Value

An object of class DNAStringSet (see the
[`Biostrings::DNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
function)

## See also

[`Biostrings::DNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)

## Author

Adrien Taudière

## Examples

``` r
dna <- physeq_or_string_to_dna(data_fungi)
dna
#> DNAStringSet object of length 1420:
#>        width seq                                            names               
#>    [1]   312 AAATGCGATAAGTAATGTGAAT...TAGGAATACCCGCTGAACTTA ASV2
#>    [2]   301 AAATGCGATAAGTAATGTGAAT...TAGGAATACCCGCTGAACTTA ASV6
#>    [3]   349 AAATGCGATAAGTAATGTGAAT...TGGGACTACCCGCTGAACTTA ASV7
#>    [4]   357 AAATGCGATAAGTAATGTGAAT...TGGGACTACCCGCTGAACTTA ASV8
#>    [5]   300 AAATGCGATAAGTAATGTGAAT...TAGGAATACCCGCTGAACTTA ASV10
#>    ...   ... ...
#> [1416]   339 AATCGCGATATGTAGTGTGATC...GGGATTACCCGCTGAACTTAA ASV1712
#> [1417]   356 GAAACGCGATATGTAATCGTGC...TAGGACTACCCGCTGAACTTA ASV1714
#> [1418]   420 AACTGCGAAACGTCATGTGACC...TAAGGATACCCGCTGAACTTA ASV1716
#> [1419]   337 GAAGTGCGATAAGCAATGCGAA...CAGGATCACCCGCTGAACTTA ASV1719
#> [1420]   260 AAACGCGAAAAGTGTTATGATG...AAGATCACCCGCTGAACTTAA ASV1726

sequences_ex <- c(
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTAATAACGAATTCATTGAATCA",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTT",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACC",
  "CGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCATGAGGTG",
  "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGG",
  "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
  "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
  "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG"
)
dna2 <- physeq_or_string_to_dna(dna_seq = sequences_ex)
dna2
#> DNAStringSet object of length 9:
#>     width seq
#> [1]    73 TACCTATGTTGCCTTGGCGGCTAAACCTACCCGG...ATGGGGCGAATTAATAACGAATTCATTGAATCA
#> [2]    70 TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTT
#> [3]    66 TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAGAGGTG
#> [4]    31 TACCTATGTTGCCTTGGCGGCTAAACCTACC
#> [5]    58 CGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCATGAGGTG
#> [6]    40 GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGG
#> [7]    40 GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG
#> [8]    40 GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG
#> [9]    40 GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG
```
