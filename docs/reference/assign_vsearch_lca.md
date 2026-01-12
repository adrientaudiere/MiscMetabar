# Assign taxonomy using LCA

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Please cite [Vsearch](https://github.com/torognes/vsearch) and
[stampa](https://github.com/frederic-mahe/stampa) if you use this
function to assign taxonomy.

1.  If top_hits_only is TRUE, the algorithm is the one of
    [stampa](https://github.com/frederic-mahe/stampa).

2.  If top_hits_only is FALSE and vote_algorithm is NULL, you need to
    carefully define `maxaccept`, `id` and `lca_cutoff` parameters. The
    algorithm is internal to vsearch using the lcaout output.

3.  If top_hits_only is FALSE and vote_algorithm is not NULL, conflict
    among the list of taxonomic assignations is resolve using the
    function
    [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md).
    The possible values for vote_algorithm are "consensus",
    "rel_majority", "abs_majority" and "unanimity". See
    [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
    for more details.

## Usage

``` r
assign_vsearch_lca(
  physeq = NULL,
  ref_fasta = NULL,
  seq2search = NULL,
  behavior = c("return_matrix", "add_to_phyloseq", "return_cmd"),
  vsearchpath = "vsearch",
  clean_pq = TRUE,
  taxa_ranks = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  nproc = 1,
  suffix = "_sintax",
  id = 0.5,
  lca_cutoff = 1,
  maxrejects = 32,
  top_hits_only = TRUE,
  maxaccepts = 0,
  keep_temporary_files = FALSE,
  verbose = TRUE,
  temporary_fasta_file = "temp.fasta",
  cmd_args = "",
  too_few = "align_start",
  vote_algorithm = NULL,
  nb_voting = NULL,
  strict = FALSE,
  nb_agree_threshold = 1,
  preference_index = NULL,
  collapse_string = "/",
  replace_collapsed_rank_by_NA = TRUE,
  simplify_taxo = TRUE,
  keep_vsearch_score = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ref_fasta:

  (required) A link to a database in vsearch format The reference
  database must contain taxonomic information in the header of each
  sequence in the form of a string starting with ";tax=" and followed by
  a comma-separated list of up to nine taxonomic identifiers. Each
  taxonomic identifier must start with an indication of the rank by one
  of the letters d (for domain) k (kingdom), p (phylum), c (class), o
  (order), f (family), g (genus), s (species), or t (strain). The letter
  is followed by a colon (:) and the name of that rank. Commas and
  semicolons are not allowed in the name of the rank. Non-ascii
  characters should be avoided in the names.

  Example:

  \\X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli,t:str.\_K-12_substr.\_MG1655

- seq2search:

  A DNAStringSet object of sequences to search for. Replace the physeq
  object.

- behavior:

  Either "return_matrix" (default), "return_cmd", or "add_to_phyloseq":

  - "return_matrix" return a list of two matrix with taxonomic value in
    the first element of the list and bootstrap value in the second one.

  - "return_cmd" return the command to run without running it.

  - "add_to_phyloseq" return a phyloseq object with amended slot
    `@taxtable`. Only available if using physeq input and not seq2search
    input.

- vsearchpath:

  (default: "vsearch") path to vsearch

- clean_pq:

  (logical, default TRUE) If set to TRUE, empty samples and empty ASV
  are discarded before clustering.

- taxa_ranks:

  A list with the name of the taxonomic rank present in ref_fasta

- nproc:

  (int, default: 1) Set to number of cpus/processors to use

- suffix:

  (character) The suffix to name the new columns. If set to "" (the
  default), the taxa_ranks algorithm is used without suffix.

- id:

  (Float \[0:1\] default 0.5). Default value is based on
  [stampa](https://github.com/frederic-mahe/stampa). See Vsearch Manual
  for parameter `--id`

- lca_cutoff:

  (int, default 1). Fraction of matching hits required for the last
  common ancestor (LCA) output. For example, a value of 0.9 imply that
  if less than 10% of assigned species are not congruent the taxonomy is
  filled. Default value is based on
  [stampa](https://github.com/frederic-mahe/stampa). See Vsearch Manual
  for parameter `--lca_cutoff`

  Text from vsearch manual : "Adjust the fraction of matching hits
  required for the last common ancestor (LCA) output with the –lcaout
  option during searches. The default value is 1.0 which requires all
  hits to match at each taxonomic rank for that rank to be included. If
  a lower cutoff value is used, e.g. 0.95, a small fraction of
  non-matching hits are allowed while that rank will still be reported.
  The argument to this option must be larger than 0.5, but not larger
  than 1.0"

- maxrejects:

  (int, default: 32) Maximum number of non-matching target sequences to
  consider before stopping the search for a given query. Default value
  is based on [stampa](https://github.com/frederic-mahe/stampa) See
  Vsearch Manual for parameter `--maxrejects`.

- top_hits_only:

  (Logical, default TRUE) Only the top hits with an equally high
  percentage of identity between the query and database sequence sets
  are written to the output. If you set top_hits_only you may need to
  set a lower `maxaccepts` and/or `lca_cutoff`. Default value is based
  on [stampa](https://github.com/frederic-mahe/stampa) See Vsearch
  Manual for parameter `--top_hits_only`

- maxaccepts:

  (int, default: 0) Default value is based on
  [stampa](https://github.com/frederic-mahe/stampa). Maximum number of
  matching target sequences to accept before stopping the search for a
  given query. See Vsearch Manual for parameter `--maxaccepts`

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files?

  - temporary_fasta_file (default "temp.fasta") : the fasta file from
    physeq or seq2search

  - "out_lca.txt" : see Vsearch Manual for parameter –lcaout

  - "userout.txt" : see Vsearch Manual for parameter –userout

- verbose:

  (logical). If TRUE, print additional information.

- temporary_fasta_file:

  Name of the temporary fasta file. Only useful with
  keep_temporary_files = TRUE.

- cmd_args:

  Additional arguments passed on to vsearch usearch_global cmd.

- too_few:

  (default value "align_start") see
  [`tidyr::separate_wider_delim()`](https://tidyr.tidyverse.org/reference/separate_wider_delim.html)

- vote_algorithm:

  (default NULL) the method to vote among "consensus", "rel_majority",
  "abs_majority" and "unanimity". See
  [`resolve_vector_ranks()`](https://adrientaudiere.github.io/MiscMetabar/reference/resolve_vector_ranks.md)
  for more details.

- nb_voting:

  (Int, default NULL). The number of taxa to keep before apply a vote to
  resolve conflict. If NULL all taxa passing the filters (min_id,
  min_bit_score, min_cover and min_e_value) are selected.

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

- keep_vsearch_score:

  (Logical, default FALSE). If TRUE, the mean and sd of id score are
  stored in the tax_table.

## Value

See param behavior

## Details

This function is mainly a wrapper of the work of others. Please cite
[vsearch](https://github.com/torognes/vsearch) and
[stampa](https://github.com/frederic-mahe/stampa)

## See also

[`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md),
[`add_new_taxonomy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/add_new_taxonomy_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
data_fungi_mini_new <- assign_vsearch_lca(data_fungi_mini,
  ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
  lca_cutoff = 0.9, behavior = "add_to_phyloseq"
)

data_fungi_mini_new2 <- assign_vsearch_lca(data_fungi_mini,
  ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
  id = 0.8, behavior = "add_to_phyloseq", top_hits_only = FALSE
)

data_fungi_mini_new3 <- assign_vsearch_lca(data_fungi_mini,
  ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
  id = 0.5, behavior = "add_to_phyloseq", top_hits_only = FALSE, vote_algorithm = "rel_majority"
)
# }
```
