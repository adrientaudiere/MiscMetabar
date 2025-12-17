# Recluster sequences of an object of class `physeq` or cluster a list of DNA sequences using vsearch software

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper of VSEARCH software.

## Usage

``` r
vsearch_clustering(
  physeq = NULL,
  dna_seq = NULL,
  nproc = 1,
  id = 0.97,
  vsearchpath = "vsearch",
  tax_adjust = 0,
  rank_propagation = FALSE,
  vsearch_cluster_method = "--cluster_size",
  vsearch_args = "--strand both",
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
  `physeq@refseq`

- nproc:

  (default: 1) Set to number of cpus/processors to use for the
  clustering

- id:

  (default: 0.97) level of identity to cluster

- vsearchpath:

  (default: "vsearch") path to vsearch

- tax_adjust:

  (Default 0) See the man page of
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  for more details. To conserved the taxonomic rank of the most abundant
  ASV, set tax_adjust to 0 (default). For the moment only tax_adjust = 0
  is robust

- rank_propagation:

  (logical, default FALSE). Do we propagate the NA value from lower
  taxonomic rank to upper rank? See the man page of
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  for more details.

- vsearch_cluster_method:

  (default: "–cluster_size) See other possible methods in the [vsearch
  manual](https://github.com/torognes/vsearch) (e.g. `--cluster_size` or
  `--cluster_smallmem`)

  - `--cluster_fast` : Clusterize the fasta sequences in filename,
    automatically sort by decreasing sequence length beforehand.

  - `--cluster_size` : Clusterize the fasta sequences in filename,
    automatically sort by decreasing sequence abundance beforehand.

  - `--cluster_smallmem` : Clusterize the fasta sequences in filename
    without automatically modifying their order beforehand. Sequence are
    expected to be sorted by decreasing sequence length, unless
    *–usersort* is used. In that case you may set `vsearch_args` to
    vsearch_args = "–strand both –usersort"

- vsearch_args:

  (default : "–strand both") a one length character element defining
  other parameters to passed on to vsearch.

- keep_temporary_files:

  (logical, default: FALSE) Do we keep temporary files ?

  - temp.fasta (refseq in fasta or dna_seq sequences)

  - cluster.fasta (centroid if method = "vsearch")

  - temp.uc (clusters if method = "vsearch")

## Value

A new object of class `physeq` or a list of cluster if dna_seq args was
used.

## Details

This function use the
[`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
function to merge taxa into clusters. By default tax_adjust = 0. See the
man page of
[`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md).

This function is mainly a wrapper of the work of others. Please cite
[vsearch](https://github.com/torognes/vsearch).

## References

VSEARCH can be downloaded from <https://github.com/torognes/vsearch>.
More information in the associated publication
<https://pubmed.ncbi.nlm.nih.gov/27781170>.

## See also

[`postcluster_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/postcluster_pq.md),
[`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/swarm_clustering.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
summary_plot_pq(data_fungi)

d_vs <- vsearch_clustering(data_fungi)
summary_plot_pq(d_vs)

# }
```
