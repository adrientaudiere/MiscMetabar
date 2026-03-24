# Recluster sequences of a phyloseq object or cluster a list of DNA sequences using MMseqs2 software

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper of the [MMseqs2](https://github.com/soedinglab/MMseqs2)
`easy-cluster` or `easy-linclust` workflow.

## Usage

``` r
mmseqs2_clustering(
  physeq = NULL,
  dna_seq = NULL,
  nproc = 1,
  id = 0.97,
  mmseqs2path = find_mmseqs2(),
  tax_adjust = 0,
  rank_propagation = FALSE,
  mmseqs2_cluster_method = c("easy-cluster", "easy-linclust"),
  coverage = 0.8,
  cov_mode = 0,
  cluster_mode = 0,
  mmseqs2_args = "",
  keep_temporary_files = FALSE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- dna_seq:

  You may directly use a character vector of DNA sequences in place of
  physeq args. When physeq is set, dna sequences take the value of
  `physeq@refseq`.

- nproc:

  (default: 1) Number of threads.

- id:

  (default: 0.97) Minimum sequence identity threshold (0–1).

- mmseqs2path:

  Path to the `mmseqs` binary (default:
  [`find_mmseqs2()`](https://adrientaudiere.github.io/MiscMetabar/reference/find_mmseqs2.md)).

- tax_adjust:

  (Default 0) See the man page of
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  for more details. To conserve the taxonomic rank of the most abundant
  taxa (ASV, OTU, ...), set `tax_adjust` to 0 (default).

- rank_propagation:

  (logical, default FALSE). Do we propagate the NA value from lower
  taxonomic rank to upper rank? See the man page of
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  for more details.

- mmseqs2_cluster_method:

  (default: `"easy-cluster"`) Either `"easy-cluster"` (cascaded
  clustering, more sensitive) or `"easy-linclust"` (linear-time
  clustering, faster for huge datasets).

- coverage:

  (numeric, default: 0.8) Alignment coverage threshold (0–1), passed to
  `-c`.

- cov_mode:

  (integer, default: 0) Coverage mode:

  - `0`: `alnRes / max(qLen, tLen)`

  - `1`: `alnRes / tLen`

  - `2`: `alnRes / qLen`

- cluster_mode:

  (integer, default: 0) Clustering algorithm:

  - `0`: greedy set cover (default)

  - `1`: connected components

  - `2`: greedy incremental

- mmseqs2_args:

  (character, default: `""`) Additional arguments passed to the MMseqs2
  clustering command.

- keep_temporary_files:

  (logical, default: FALSE) Keep intermediate files for debugging?

## Value

A new object of class `physeq` or a data.frame of cluster membership if
`dna_seq` was used.

## Details

This function is mainly a wrapper of the work of others. Please cite
[MMseqs2](https://github.com/soedinglab/MMseqs2).

## References

MMseqs2 can be downloaded from <https://github.com/soedinglab/MMseqs2>.
More information in the associated publication
[doi:10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988) .

## See also

[`postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md),
[`vsearch_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/vsearch_clustering.md),
[`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/swarm_clustering.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
d_mm <- mmseqs2_clustering(data_fungi_mini)
# }
```
