# Assign taxonomy using blastn algorithm and the blast software

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the blast software.

## Usage

``` r
assign_blastn(
  physeq,
  ref_fasta = NULL,
  database = NULL,
  blastpath = NULL,
  behavior = c("return_matrix", "add_to_phyloseq"),
  method_algo = c("vote", "top-hit"),
  suffix = "_blastn",
  min_id = 95,
  min_bit_score = 50,
  min_cover = 95,
  min_e_value = 1e-30,
  nb_voting = NULL,
  column_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  vote_algorithm = c("consensus", "rel_majority", "abs_majority", "unanimity"),
  strict = FALSE,
  nb_agree_threshold = 1,
  preference_index = NULL,
  collapse_string = "/",
  replace_collapsed_rank_by_NA = TRUE,
  simplify_taxo = TRUE,
  keep_blast_metrics = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ref_fasta:

  Either a DNAStringSet object or a path to a fasta file to make the
  blast database. It must be in sintax format. See
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md).

- database:

  path to a blast database. Only used if ref_fasta is not set.

- blastpath:

  path to blast program.

- behavior:

  Either "return_matrix" (default), or "add_to_phyloseq":

  - "return_matrix" return a list of two matrix with taxonomic value in
    the first element of the list and bootstrap value in the second one.

  - "add_to_phyloseq" return a phyloseq object with amended slot
    `@taxtable`. Only available if using physeq input and not seq2search
    input.

- method_algo:

  (One of "vote" or "top-hit"). If top-hit, only the better match is
  used to assign taxonomy. If vote, the algorithm takes all (or
  `nb_voting` if `nb_voting` is not null) select assignation and resolve
  the conflict using the function
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md).

- suffix:

  (character) The suffix to name the new columns. If set to "" (the
  default), the taxa_ranks algorithm is used without suffix.

- min_id:

  (default: 95) the identity percent to take into account a references
  taxa

- min_bit_score:

  (default: 50) the minimum bit score to take into account a references
  taxa

- min_cover:

  (default: 50) cut of in query cover (%) to keep result

- min_e_value:

  (default: 1e-30) cut of in e-value (%) to keep result The BLAST
  E-value is the number of expected hits of similar quality (score) that
  could be found just by chance.

- nb_voting:

  (Int, default NULL). The number of taxa to keep before apply a vote to
  resolve conflict. If NULL all taxa passing the filters (min_id,
  min_bit_score, min_cover and min_e_value) are selected.

- column_names:

  A vector of names for taxonomic ranks. Must correspond to names in the
  ref_fasta files.

- vote_algorithm:

  the method to vote among "consensus", "rel_majority", "abs_majority"
  and "unanimity". See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  for more details.

- strict:

  (Logical, default FALSE). See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  for more details.

- nb_agree_threshold:

  See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  for more details.

- preference_index:

  See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  for more details.

- collapse_string:

  See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  for more details.

- replace_collapsed_rank_by_NA:

  (Logical, default TRUE) See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  for more details.

- simplify_taxo:

  (logical default TRUE). Do we apply the function
  [`simplify_taxo()`](https://adrientaudiere.github.io/MiscMetabar/reference/simplify_taxo.md)
  to the phyloseq object?

- keep_blast_metrics:

  (Logical, default FALSE). If TRUE, the blast metrics ("Query seq.
  length", "Taxa seq. length", "Alignment length", "% id. match",
  "e-value", "bit score" and "Query cover") are stored in the tax_table.

- ...:

  Additional arguments passed on to
  [`blast_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/blast_pq.md)

## Value

- If behavior == "return_matrix" :

  - If method_algo = "top-hit" a matrix of taxonomic assignation

  - If method_algo = "vote", a list of two matrix, the first is the raw
    taxonomic assignation (before vote). The second one is the taxonomic
    assignation in which conflicts are resolved using vote.

- If behavior == "add_to_phyloseq", return a new phyloseq object

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
ref_fasta <- Biostrings::readDNAStringSet(system.file("extdata",
  "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar", mustWork = TRUE
))

# assign_blastn(data_fungi_mini, ref_fasta = ref_fasta) # error because not
# enough sequences in db so none blast query passed the filters.
# So we used low score filter hereafter.

mat <- assign_blastn(data_fungi_mini,
  ref_fasta = ref_fasta,
  method_algo = "top-hit", min_id = 70, min_e_value = 1e-3, min_cover = 50,
  min_bit_score = 20
)
head(mat)

assign_blastn(data_fungi_mini,
  ref_fasta = ref_fasta, method_algo = "vote",
  vote_algorithm = "rel_majority", min_id = 90, min_cover = 50,
  behavior = "add_to_phyloseq"
)@tax_table

assign_blastn(data_fungi_mini,
  ref_fasta = ref_fasta, method_algo = "vote",
  vote_algorithm = "consensus", replace_collapsed_rank_by_NA = FALSE,
  min_id = 90, min_cover = 50, behavior = "add_to_phyloseq"
)@tax_table
} # }
```
