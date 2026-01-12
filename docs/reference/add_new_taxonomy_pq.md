# Add new taxonomic rank to a phyloseq object.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

One of main use of this function is to add taxonomic assignment from a
new database.

## Usage

``` r
add_new_taxonomy_pq(
  physeq,
  ref_fasta,
  suffix = NULL,
  method = c("dada2", "sintax", "lca", "idtaxa", "blastn", "dada2_2steps"),
  trainingSet = NULL,
  min_bootstrap = NULL,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- ref_fasta:

  (required) A link to a database. passed on to
  [`dada2::assignTaxonomy`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html).

- suffix:

  (character) The suffix to name the new columns. If set to NULL (the
  default), the basename of the file reFasta is used with the name of
  the method. Set suffix to "" in order to remove any suffix.

- method:

  (required, default "dada2") :

  - "dada2":
    [`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html)

  - "dada2_2step":
    [`assign_dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_dada2.md)

  - "sintax": see
    [`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md)

  - "lca": see
    [`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md)

  - "idtaxa": see
    [`assign_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_idtaxa.md)

  - "blastn": see
    [`assign_blastn()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_blastn.md)

- trainingSet:

  see
  [`assign_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_idtaxa.md).
  Only used if method = "idtaxa". Note that if trainingSet is not NULL,
  the ref_fasta is overwrite by the trainingSet parameter. To customize
  learning parameters of the idtaxa algorithm you must use trainingSet
  computed by the function
  [`learn_idtaxa()`](https://adrientaudiere.github.io/MiscMetabar/reference/learn_idtaxa.md).

- min_bootstrap:

  (Float \[0:1\])

  Minimum bootstrap value to inform taxonomy. For each bootstrap below
  the min_bootstrap value, the taxonomy information is set to NA.

  Correspond to parameters :

  - dada2 & dada2_2step: `minBoot`, default value = 0.5

  - sintax: `min_bootstrap`, default value = 0.5

  - lca: `id`, default value = 0.5. Note in that case, the bootstrap
    value is different. See the id parameter in
    [`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md)

  - idtaxa: `threshold`, default value = 0.6

  - blastn: This method do not take different bootstrap value. You may
    use method="vote" with different `vote_algorithm` as well as
    different filters parameters (min_id, min_bit_score, min_cover and
    min_e_value)

- ...:

  Additional arguments passed on to the taxonomic assignation method.

## Value

A new
[`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with a larger slot tax_table"

## See also

[`dada2::assignTaxonomy()`](https://rdrr.io/pkg/dada2/man/assignTaxonomy.html),
[`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md),
[`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md),
[`assign_sintax()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_sintax.md),
[`assign_blastn()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_blastn.md),
[`assign_dada2()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_dada2.md)

## Author

Adrien Taudière

## Examples

``` r
if (FALSE) { # \dontrun{
ref_fasta <- system.file("extdata",
  "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar", mustWork = TRUE
)
add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "dada2")
add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "dada2_2steps")
add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "lca")
add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "idtaxa")

# blastn doesn't work with fasta.gz format
ref_fasta <- system.file("extdata",
  "100_sp_UNITE_sh_general_release_dynamic_sintax.fasta",
  package = "MiscMetabar", mustWork = TRUE
)

dp <- add_new_taxonomy_pq(data_fungi_mini, ref_fasta,
  method = "blastn", min_id = 80, min_cover = 50, min_bit_score = 20,
  min_e_value = 1e-20
)
dp_tophit <- add_new_taxonomy_pq(data_fungi_mini, ref_fasta,
  method = "blastn", min_id = 80, min_cover = 50, min_bit_score = 20,
  min_e_value = 1e-20, method_algo = "top_hit"
)
} # }
```
