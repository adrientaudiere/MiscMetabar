# Recluster sequences of an object of class `physeq` or a list of DNA sequences

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function use the `merge_taxa_vec` function to merge taxa into
clusters.

## Usage

``` r
postcluster_pq(
  physeq = NULL,
  dna_seq = NULL,
  nproc = 1,
  method = "clusterize",
  id = 0.97,
  vsearchpath = "vsearch",
  tax_adjust = 0,
  rank_propagation = FALSE,
  vsearch_cluster_method = "--cluster_size",
  vsearch_args = "--strand both",
  keep_temporary_files = FALSE,
  swarmpath = "swarm",
  d = 1,
  swarm_args = "--fastidious",
  method_clusterize = "overlap",
  ...
)

asv2otu(
  physeq = NULL,
  dna_seq = NULL,
  nproc = 1,
  method = "clusterize",
  id = 0.97,
  vsearchpath = "vsearch",
  tax_adjust = 0,
  rank_propagation = FALSE,
  vsearch_cluster_method = "--cluster_size",
  vsearch_args = "--strand both",
  keep_temporary_files = FALSE,
  swarmpath = "swarm",
  d = 1,
  swarm_args = "--fastidious",
  method_clusterize = "overlap",
  ...
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

- method:

  (default: clusterize) Set the clustering method.

  - `clusterize` use the
    [`DECIPHER::Clusterize()`](https://rdrr.io/pkg/DECIPHER/man/Clusterize.html)
    fonction,

  - `vsearch` use the vsearch software
    (https://github.com/torognes/vsearch) with arguments
    `--cluster_size` by default (see args `vsearch_cluster_method`) and
    `-strand both` (see args `vsearch_args`)

  - `swarm` use the swarm

- id:

  (default: 0.97) level of identity to cluster

- vsearchpath:

  (default: vsearch) path to vsearch

- tax_adjust:

  (Default 0) See the man page of
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  for more details. To conserved the taxonomic rank of the most abundant
  taxa (ASV, OTU,...), set tax_adjust to 0 (default). For the moment
  only tax_adjust = 0 is robust

- rank_propagation:

  (logical, default FALSE). Do we propagate the NA value from lower
  taxonomic rank to upper rank? See the man page of
  [`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md)
  for more details.

- vsearch_cluster_method:

  (default: "–cluster_size) See other possible methods in the [vsearch
  manual](https://github.com/torognes/vsearch/) (e.g. `--cluster_size`
  or `--cluster_smallmem`)

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

  (logical, default: FALSE) Do we keep temporary files

  - temp.fasta (refseq in fasta or dna_seq sequences)

  - cluster.fasta (centroid if method = "vsearch")

  - temp.uc (clusters if method = "vsearch")

- swarmpath:

  (default: swarm) path to swarm

- d:

  (default: 1) maximum number of differences allowed between two
  amplicons, meaning that two amplicons will be grouped if they have `d`
  (or less) differences

- swarm_args:

  (default : "–fastidious") a one length character element defining
  other parameters to passed on to swarm See other possible methods in
  the [SWARM pdf
  manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)

- method_clusterize:

  (default "overlap") the method for the
  [`DECIPHER::Clusterize()`](https://rdrr.io/pkg/DECIPHER/man/Clusterize.html)
  method

- ...:

  Additional arguments passed on to
  [`DECIPHER::Clusterize()`](https://rdrr.io/pkg/DECIPHER/man/Clusterize.html)

## Value

A new object of class `physeq` or a list of cluster if dna_seq args was
used.

## Details

This function use the `merge_taxa_vec` function to merge taxa into
clusters. By default tax_adjust = 0. See the man page of
[`merge_taxa_vec()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_taxa_vec.md).

## References

VSEARCH can be downloaded from <https://github.com/torognes/vsearch>.
More information in the associated publication
<https://pubmed.ncbi.nlm.nih.gov/27781170>.

## See also

[`vsearch_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/vsearch_clustering.md)
and
[`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/reference/swarm_clustering.md)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("DECIPHER")) {
  postcluster_pq(data_fungi_mini)
}
#> Loading required namespace: DECIPHER
#> Partitioning sequences by 3-mer similarity:
#> ================================================================================
#> 
#> Time difference of 0.03 secs
#> 
#> Sorting by relatedness within 11 groups:
#> 
iteration 1 of up to 17 (100.0% stability) 
#> 
#> Time difference of 0.01 secs
#> 
#> Clustering sequences by 9-mer similarity:
#> ================================================================================
#> 
#> Time difference of 0.08 secs
#> 
#> Clusters via relatedness sorting: 100% (0% exclusively)
#> Clusters via rare 3-mers: 100% (0% exclusively)
#> Estimated clustering effectiveness: 100%
#> 
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 32 taxa and 137 samples ]
#> sample_data() Sample Data:       [ 137 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 32 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 32 reference sequences ]
# \donttest{
if (requireNamespace("DECIPHER")) {
  postcluster_pq(data_fungi_mini, method_clusterize = "longest")

  if (MiscMetabar::is_swarm_installed()) {
    d_swarm <- postcluster_pq(data_fungi_mini, method = "swarm")
  }
  if (MiscMetabar::is_vsearch_installed()) {
    d_vs <- postcluster_pq(data_fungi_mini, method = "vsearch")
  }
}
#> Partitioning sequences by 3-mer similarity:
#> ================================================================================
#> 
#> Time difference of 0.02 secs
#> 
#> Sorting by relatedness within 11 groups:
#> 
iteration 1 of up to 17 (100.0% stability) 
#> 
#> Time difference of 0.01 secs
#> 
#> Clustering sequences by 9-mer similarity:
#> ================================================================================
#> 
#> Time difference of 0.08 secs
#> 
#> Clusters via relatedness sorting: 100% (0% exclusively)
#> Clusters via rare 3-mers: 100% (0% exclusively)
#> Estimated clustering effectiveness: 100%
#> 
# }
```
