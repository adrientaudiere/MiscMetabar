# Cluster raw reads into OTUs (vsearch or SWARM), no DADA2 required

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Cluster marker-gene reads into similarity-threshold OTUs and assemble
them into a
[`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with an `otu_table` of per-sample abundances and the
representative (seed/centroid) sequence of each OTU in the `refseq`
slot: reads are dereplicated per sample with
[vsearch](https://github.com/torognes/vsearch), pooled, dereplicated
again, and clustered with vsearch (`--cluster_size`, default) or
[SWARM](https://github.com/torognes/swarm) (`d > 1` recommended).

This complements
[`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/postcluster_pq.md),
which post-clusters ASVs already stored in a phyloseq object:
`cluster_reads()` starts one step earlier, directly from sequence files
or a vector of DNA sequences. Use
[`denoised_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/denoised_reads.md)
instead of this function to denoise reads into ASVs/zOTUs (dada2
pipeline, SWARM `d = 1` or vsearch UNOISE).

## Usage

``` r
cluster_reads(
  path_to_fastx = NULL,
  dna_seq = NULL,
  sample_name = "sample_1",
  method = c("vsearch", "swarm"),
  id = 0.97,
  d = 1,
  fastidious = TRUE,
  minsize = 1,
  nproc = 1,
  vsearchpath = find_vsearch(),
  swarmpath = "swarm",
  vsearch_args = "--strand both",
  swarm_args = "",
  keep_temporary_files = FALSE
)
```

## Arguments

- path_to_fastx:

  (character) Path(s) to fasta or fastq files (optionally gzipped). Each
  file is treated as one sample; the sample name is the basename
  stripped of its extension. Set either `path_to_fastx` or `dna_seq`,
  not both.

- dna_seq:

  A named character vector or
  [`Biostrings::DNAStringSet-class()`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  of DNA sequences, treated as one sample. Duplicated sequences are
  counted; alternatively, already dereplicated sequences may carry their
  abundance as a `;size=N` suffix in their names. Set either
  `path_to_fastx` or `dna_seq`, not both.

- sample_name:

  (default: `"sample_1"`) Sample name used when `dna_seq` is set.

- method:

  (default: `"vsearch"`) Clustering tool: `"vsearch"` or `"swarm"`.

- id:

  (default: 0.97) Identity threshold for vsearch clustering.

- d:

  (default: 1) SWARM resolution. `d = 1` is a denoising resolution and
  triggers a warning pointing to
  [`denoised_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/denoised_reads.md);
  use `d > 1` for OTU-like clustering.

- fastidious:

  (logical, default: TRUE) Use SWARM fastidious mode. Ignored when
  `d != 1`.

- minsize:

  (default: 1) Discard dereplicated sequences with an abundance strictly
  lower than `minsize` before clustering. The common pre-clustering
  filter `minsize = 2` removes singletons.

- nproc:

  (default: 1) Number of threads for the clustering tool.

- vsearchpath:

  (default:
  [`find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_vsearch.md))
  Path to the vsearch binary.

- swarmpath:

  (default: `"swarm"`) Path to the swarm binary.

- vsearch_args:

  (default: `"--strand both"`) Extra arguments passed to vsearch
  clustering.

- swarm_args:

  (default: `""`) Extra arguments passed to swarm.

- keep_temporary_files:

  (logical, default: FALSE) Keep the dereplicated fasta files, the
  pooled fasta, the cluster table and the representative-sequences fasta
  in a temporary directory.

## Value

A
[`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with one row (taxon) per OTU named `OTU_1`, `OTU_2`, ... ordered
by decreasing total abundance, one column (sample) per input file, the
representative sequence of each OTU in the `refseq` slot, and (for file
input) the source file paths in the `sample_data` slot.

## Details

vsearch is required (per-sample and pool-level dereplication, and
clustering for `method = "vsearch"`); swarm is additionally required for
`method = "swarm"`. Use
[`is_vsearch_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_vsearch_installed.md)
and
[`is_swarm_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_swarm_installed.md)
to check availability. This function is mainly a wrapper of the work of
others: please cite [vsearch](https://github.com/torognes/vsearch)
(<https://pubmed.ncbi.nlm.nih.gov/27781170>) and/or
[SWARM](https://github.com/torognes/swarm)
([doi:10.1093/bioinformatics/btab493](https://doi.org/10.1093/bioinformatics/btab493)
).

## See also

[`denoised_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/denoised_reads.md),
[`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/postcluster_pq.md),
[`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/swarm_clustering.md),
[`vsearch_clustering()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/vsearch_clustering.md),
[`find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_vsearch.md),
[`install_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/install_vsearch.md)

## Author

Adrien Taudière

## Examples

``` r
f1 <- system.file("extdata", "ex_R1_001.fasta", package = "MiscMetabar")
f2 <- system.file("extdata", "ex_R1_002.fasta", package = "MiscMetabar")
(pq_otu <- cluster_reads(path_to_fastx = c(f1, f2)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 50 taxa and 2 samples ]
#> sample_data() Sample Data:       [ 2 samples by 1 sample variables ]
#> refseq()      DNAStringSet:      [ 50 reference sequences ]

cluster_reads(path_to_fastx = c(f1, f2), method = "swarm", d = 3)
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 50 taxa and 2 samples ]
#> sample_data() Sample Data:       [ 2 samples by 1 sample variables ]
#> refseq()      DNAStringSet:      [ 50 reference sequences ]
```
