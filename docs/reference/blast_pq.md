# Blast all sequence of `refseq` slot of a [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object against a custom database.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the blast software.

## Usage

``` r
blast_pq(
  physeq,
  fasta_for_db = NULL,
  database = NULL,
  blastpath = NULL,
  id_cut = 90,
  bit_score_cut = 50,
  min_cover_cut = 50,
  e_value_cut = 1e-30,
  unique_per_seq = FALSE,
  score_filter = TRUE,
  nproc = 1,
  args_makedb = NULL,
  args_blastn = NULL,
  keep_temporary_files = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fasta_for_db:

  Either a DNAStringSet object or a path to a fasta file to make the
  blast database.

- database:

  path to a blast database

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

- nproc:

  (default: 1) Set to number of cpus/processors to use for blast (args
  -num_threads for blastn command)

- args_makedb:

  Additional arguments passed on to makeblastdb command

- args_blastn:

  Additional arguments passed on to blastn command

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files

  - db.fasta (refseq transformed into a database)

  - dbase list of files (output of blastn)

  - blast_result.txt the summary result of blastn using
    `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`

  - temp.fasta if fasta_for_db was a DNAStringSet object.

## Value

a blast table

## See also

[`blast_to_phyloseq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_to_phyloseq.md)
to use `refseq` slot as a database
