# Search for a list of sequence in a fasta file against physeq reference sequences using [vsearch](https://github.com/torognes/vsearch)

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use of VSEARCH software.

## Usage

``` r
vs_search_global(
  physeq,
  path_to_fasta = NULL,
  seq2search = NULL,
  vsearchpath = "vsearch",
  id = 0.8,
  iddef = 0,
  keep_temporary_files = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- path_to_fasta:

  (required if seq2search is NULL) a path to fasta file if seq2search is
  est to NULL.

- seq2search:

  (required if path_to_fasta is NULL) Either (i) a DNAStringSet object
  or (ii) a character vector that will be convert to DNAStringSet using
  [`Biostrings::DNAStringSet()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)

- vsearchpath:

  (default: "vsearch") path to vsearch

- id:

  (default: 0.8) id for the option `--usearch_global` of the vsearch
  software

- iddef:

  (default: 0) iddef for the option `--usearch_global` of the vsearch
  software

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files

  - temp.fasta (refseq in fasta)

  - cluster.fasta (centroid)

  - temp.uc (clusters)

## Value

A dataframe with uc results (invisible)

## Details

This function is mainly a wrapper of the work of others. Please cite
[vsearch](https://github.com/torognes/vsearch).

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("seqinr")) {
  file_dna <- tempfile("dna.fa")
  seqinr::write.fasta("GCCCATTAGTATTCTAGTGGGCATGCCTGTTCGAGCGTCATTTTCAACC",
    file = file_dna, names = "seq1"
  )

  res <- vs_search_global(data_fungi, path_to_fasta = file_dna)
  unlink(file_dna)

  res[res$identity != "*", ]

  clean_pq(subset_taxa(data_fungi, res$identity != "*"))
}
#> Error in eval(e, x, parent.frame()): object 'res' not found
# }
```
