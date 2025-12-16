# Search for exact matching of sequences

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Search for exact matching of sequences using complement, reverse and
reverse-complement. It is useful to check for primers issues after
cutadapt step.

## Usage

``` r
search_exact_seq_pq(physeq, seq2search)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- seq2search:

  A DNAStringSet object of sequences to search for.

## Value

A list of data-frames for each input sequences with the name, the
sequences and the number of occurrences of the original sequence, the
complement sequence, the reverse sequence and the reverse-complement
sequence.

## Author

Adrien Taudière

## Examples

``` r
data("data_fungi")
search_primers <- search_exact_seq_pq(data_fungi,
  seq2search = Biostrings::DNAStringSet(c("TTGAACGCACATTGCGCC", "ATCCCTACCTGATCCGAG"))
)
```
