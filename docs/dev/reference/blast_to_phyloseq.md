# Blast some sequence against `refseq` slot of a [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the blast software.

## Usage

``` r
blast_to_phyloseq(
  physeq,
  seq2search,
  blastpath = NULL,
  id_cut = 90,
  bit_score_cut = 50,
  min_cover_cut = 50,
  e_value_cut = 1e-30,
  unique_per_seq = FALSE,
  score_filter = TRUE,
  list_no_output_query = FALSE,
  args_makedb = NULL,
  args_blastn = NULL,
  nproc = 1,
  keep_temporary_files = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

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

- args_makedb:

  Additional arguments passed on to makeblastdb command

- args_blastn:

  Additional arguments passed on to blastn command

- nproc:

  (default: 1) Set to number of cpus/processors to use for blast (args
  -num_threads for blastn command)

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files

  - db.fasta (refseq transformed into a database)

  - dbase list of files (output of blastn)

  - blast_result.txt the summary result of blastn using
    `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`

## Value

the blast table

## See also

[`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/blast_pq.md)
to use `refseq` slot as query sequences against un custom database.

## Examples

``` r
blast_to_phyloseq(data_fungi,
  seq2search = system.file("extdata", "ex.fasta",
    package = "MiscMetabar", mustWork = TRUE
  )
)
#>           Query name Query seq. length Taxa name Taxa seq. length
#> 1     Eutypa_spinosa               467      ASV2              312
#> 501  Eutypa_spinosa2               464      ASV2              312
#> 1001 Eutypa_spinosa3               465      ASV2              312
#> 2     Eutypa_spinosa               467   ASV1284              312
#> 3     Eutypa_spinosa               467    ASV375              312
#> 502  Eutypa_spinosa2               464   ASV1284              312
#> 503  Eutypa_spinosa2               464    ASV375              312
#> 1002 Eutypa_spinosa3               465   ASV1284              312
#> 1003 Eutypa_spinosa3               465    ASV375              312
#> 4     Eutypa_spinosa               467    ASV576              313
#> 504  Eutypa_spinosa2               464    ASV576              313
#> 1004 Eutypa_spinosa3               465    ASV576              313
#> 5     Eutypa_spinosa               467    ASV541              312
#> 505  Eutypa_spinosa2               464    ASV541              312
#> 1005 Eutypa_spinosa3               465    ASV541              312
#> 6     Eutypa_spinosa               467   ASV1323              312
#> 7     Eutypa_spinosa               467   ASV1082              312
#> 8     Eutypa_spinosa               467   ASV1027              312
#> 506  Eutypa_spinosa2               464   ASV1323              312
#> 507  Eutypa_spinosa2               464   ASV1082              312
#> 508  Eutypa_spinosa2               464   ASV1027              312
#> 1006 Eutypa_spinosa3               465   ASV1323              312
#> 1007 Eutypa_spinosa3               465   ASV1082              312
#> 1008 Eutypa_spinosa3               465   ASV1027              312
#> 9     Eutypa_spinosa               467   ASV1616              313
#> 10    Eutypa_spinosa               467    ASV634              312
#> 509  Eutypa_spinosa2               464   ASV1616              313
#> 510  Eutypa_spinosa2               464    ASV634              312
#> 1009 Eutypa_spinosa3               465   ASV1616              313
#> 1010 Eutypa_spinosa3               465    ASV634              312
#> 11    Eutypa_spinosa               467    ASV962              313
#> 511  Eutypa_spinosa2               464    ASV962              313
#> 1011 Eutypa_spinosa3               465    ASV962              313
#>      Alignment length % id. match  e-value bit score Query cover
#> 1                 234      93.162 1.26e-95       342          50
#> 501               234      93.162 1.25e-95       342          50
#> 1001              234      93.162 1.26e-95       342          50
#> 2                 234      92.308 2.73e-92       331          50
#> 3                 234      92.308 2.73e-92       331          50
#> 502               234      92.308 2.71e-92       331          50
#> 503               234      92.308 2.71e-92       331          50
#> 1002              234      92.308 2.72e-92       331          50
#> 1003              234      92.308 2.72e-92       331          50
#> 4                 235      91.915 3.53e-91       327          50
#> 504               235      91.915 3.51e-91       327          50
#> 1004              235      91.915 3.52e-91       327          50
#> 5                 234      91.880 1.27e-90       326          50
#> 505               234      91.880 1.26e-90       326          50
#> 1005              234      91.880 1.26e-90       326          50
#> 6                 234      91.453 5.91e-89       320          50
#> 7                 235      91.489 5.91e-89       320          50
#> 8                 234      91.453 5.91e-89       320          50
#> 506               234      91.453 5.87e-89       320          50
#> 507               235      91.489 5.87e-89       320          50
#> 508               234      91.453 5.87e-89       320          50
#> 1006              234      91.453 5.88e-89       320          50
#> 1007              235      91.489 5.88e-89       320          50
#> 1008              234      91.453 5.88e-89       320          50
#> 9                 234      91.026 2.75e-87       315          50
#> 10                234      91.026 2.75e-87       315          50
#> 509               234      91.026 2.73e-87       315          50
#> 510               234      91.026 2.73e-87       315          50
#> 1009              234      91.026 2.74e-87       315          50
#> 1010              234      91.026 2.74e-87       315          50
#> 11                235      90.213 1.66e-84       305          50
#> 511               235      90.213 1.64e-84       305          50
#> 1011              235      90.213 1.65e-84       305          50
```
