# Assign taxonomy using MMseqs2

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Use the [MMseqs2](https://github.com/soedinglab/MMseqs2) software to
assign taxonomy to sequences.

The preferred usage is to provide a reference FASTA file in SINTAX
format via `ref_fasta`. The function builds a temporary MMseqs2 taxonomy
database from the SINTAX headers and then runs `mmseqs easy-taxonomy`
with the requested `--lca-mode`, giving the same LCA behaviour as the
`database` path.

Alternatively, a pre-built MMseqs2 database with NCBI taxonomy can be
passed via the `database` parameter (created via `mmseqs createdb` +
`mmseqs createtaxdb`, or downloaded with `mmseqs databases`). In this
case, the MMseqs2 native `easy-taxonomy` LCA workflow is used. See the
[MMseqs2 wiki](https://github.com/soedinglab/MMseqs2/wiki) for details.

## Usage

``` r
assign_mmseqs2(
  physeq = NULL,
  ref_fasta = NULL,
  database = NULL,
  seq2search = NULL,
  mmseqs2path = find_mmseqs2(),
  behavior = c("return_matrix", "add_to_phyloseq"),
  suffix = "_mmseqs2",
  lca_mode = 3,
  lca_ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species"),
  column_names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  search_type = 3,
  sensitivity = NULL,
  min_seq_id = NULL,
  e_value = NULL,
  max_accept = 5,
  nproc = 1,
  clean_pq = TRUE,
  simplify_taxo = TRUE,
  keep_temporary_files = FALSE,
  verbose = FALSE,
  cmd_args = ""
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ref_fasta:

  Either a
  [Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object or a path to a FASTA file in SINTAX format (taxonomy in headers
  after `;tax=`). Only used if `database` is not set. See
  [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md)
  for the SINTAX format specification.

- database:

  (optional) Path to a pre-built MMseqs2 database with NCBI taxonomy
  information. Only used if `ref_fasta` is not set.

- seq2search:

  (optional) A
  [Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object. Use instead of `physeq` to search arbitrary sequences. Cannot
  be used together with `physeq`.

- mmseqs2path:

  Path to the `mmseqs` binary (default:
  [`find_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_mmseqs2.md)).

- behavior:

  Either `"return_matrix"` (default) or `"add_to_phyloseq"`:

  - `"return_matrix"`: return a data frame with taxonomic assignments.

  - `"add_to_phyloseq"`: return a phyloseq object with the taxonomy
    appended to the `tax_table` slot.

- suffix:

  (character) Suffix appended to new taxonomy column names (default:
  `"_mmseqs2"`).

- lca_mode:

  (integer) The LCA mode used by MMseqs2:

  - `1`: single search LCA

  - `3` (default): approximate 2bLCA (fast, recommended)

  - `4`: top-hit LCA (all equal-scoring top hits)

- lca_ranks:

  Character vector of NCBI taxonomy rank names passed to `--lca-ranks`
  (default:
  `c("superkingdom", "phylum", "class", "order", "family", "genus", "species")`).

- column_names:

  Character vector of output column names, must be the same length as
  `lca_ranks` (default:
  `c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")`).

- search_type:

  (integer) MMseqs2 search type:

  - `0`: auto-detect

  - `2`: translated nucleotide

  - `3` (default): nucleotide

- sensitivity:

  (numeric, optional) Search sensitivity (`-s` parameter). Higher values
  are slower but more sensitive (range 1–7). If `NULL`, MMseqs2 uses its
  default.

- min_seq_id:

  (numeric, optional) Minimum sequence identity (0–1). If `NULL`,
  MMseqs2 uses its default.

- e_value:

  (numeric, optional) Maximum E-value threshold (`-e`). If `NULL`,
  MMseqs2 uses its default.

- max_accept:

  (integer, optional) Maximum number of hits accepted per query
  (`--max-accept`). Useful with `lca_mode = 1` or `4` to widen the hit
  set used for LCA (default: `5`).

- nproc:

  (integer) Number of threads (default: 1).

- clean_pq:

  (logical) Clean the phyloseq object before searching? (default:
  `TRUE`).

- simplify_taxo:

  (logical) Apply
  [`simplify_taxo()`](https://adrientaudiere.github.io/MiscMetabar/reference/simplify_taxo.md)
  to the result? Only used when `behavior = "add_to_phyloseq"` (default:
  `TRUE`).

- keep_temporary_files:

  (logical) Keep intermediate files for debugging? (default: `FALSE`).

- verbose:

  (logical) Print progress messages? (default: `FALSE`).

- cmd_args:

  (character) Additional arguments appended to the MMseqs2 command.

## Value

- If `behavior == "return_matrix"`: a
  [tibble](https://tibble.tidyverse.org/reference/tibble.html) with
  columns `taxa_names` and one column per rank.

- If `behavior == "add_to_phyloseq"`: a new phyloseq object with amended
  `tax_table`.

## Details

This function is mainly a wrapper of the work of others. Please cite
[MMseqs2](https://github.com/soedinglab/MMseqs2): Mirdita M, Steinegger
M, Breitwieser F, Soding J, Levy Karin E: Fast and sensitive taxonomic
assignment to metagenomic contigs. Bioinformatics (2021).

## See also

[`assign_blastn()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_blastn.md),
[`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md),
[`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
ref_fasta <- Biostrings::readDNAStringSet(system.file("extdata",
  "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar", mustWork = TRUE
))

# Preferred usage: provide a SINTAX-formatted FASTA file.
# The function searches with easy-search and parses SINTAX headers.
res <- assign_mmseqs2(data_fungi_mini, ref_fasta = ref_fasta)
head(res)

# Add taxonomy to phyloseq:
physeq_new <- assign_mmseqs2(
  data_fungi_mini,
  ref_fasta = ref_fasta,
  behavior = "add_to_phyloseq"
)
} # }
```
