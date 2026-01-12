# Blast some sequence against sequences from of a [`derep-class`](https://rdrr.io/pkg/dada2/man/derep-class.html) object.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the blast software.

## Usage

``` r
blast_to_derep(
  derep,
  seq2search,
  blastpath = NULL,
  id_cut = 90,
  bit_score_cut = 50,
  min_cover_cut = 50,
  e_value_cut = 1e-30,
  unique_per_seq = FALSE,
  score_filter = FALSE,
  list_no_output_query = FALSE,
  min_length_seq = 200,
  args_makedb = NULL,
  args_blastn = NULL,
  nproc = 1,
  keep_temporary_files = FALSE
)
```

## Arguments

- derep:

  The result of
  [`dada2::derepFastq()`](https://rdrr.io/pkg/dada2/man/derepFastq.html).
  A list of `derep-class` object.

- seq2search:

  (required) path to a fasta file defining the sequences you want to
  blast against the taxa (ASV, OTU) sequences from the physeq object.

- blastpath:

  path to blast program

- id_cut:

  (default: 90) cut of in identity percent to keep result

- bit_score_cut:

  (default: 50) cut of in bit score to keep result The higher the
  bit-score, the better the sequence similarity. The bit-score is the
  requires size of a sequence database in which the current match could
  be found just by chance. The bit-score is a log2 scaled and normalized
  raw-score. Each increase by one doubles the required database size
  (2bit-score).

- min_cover_cut:

  (default: 50) cut of in query cover (%) to keep result

- e_value_cut:

  (default: 1e-30) cut of in e-value (%) to keep result The BLAST
  E-value is the number of expected hits of similar quality (score) that
  could be found just by chance.

- unique_per_seq:

  (logical, default FALSE) if TRUE only return the better match (higher
  **bit score**) for each sequence

- score_filter:

  (logical, default TRUE) does results are filter by score? If FALSE,
  `id_cut`,`bit_score_cut`, `e_value_cut` and `min_cover_cut` are
  ignored

- list_no_output_query:

  (logical) does the result table include query sequences for which
  `blastn` does not find any correspondence?

- min_length_seq:

  (default: 200) Removed sequences with less than `min_length_seq` from
  derep before blast. Set to 0 to discard filtering sequences by length.

- args_makedb:

  Additional arguments passed on to makeblastdb command

- args_blastn:

  Additional arguments passed on to blastn command

- nproc:

  (default: 1) Set to number of cpus/processors to use for blast (args
  -num_threads for blastn command)

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files :

  - db.fasta (refseq transformed into a database)

  - dbase list of files (output of blastn)

  - blast_result.txt the summary result of blastn using
    `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`

## Value

A blast table

## See also

[`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md)
to use `refseq` slot as query sequences against un custom database and
[`blast_to_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_phyloseq.md)
to use `refseq` slot as a database

## Author

Adrien Taudière
