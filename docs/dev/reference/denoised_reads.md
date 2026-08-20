# Denoise raw reads into ASVs / zOTUs (dada2, SWARM or vsearch UNOISE)

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Denoise marker-gene reads into amplicon sequence variants (ASVs, named
`ASV_1`, `ASV_2`, ...) and assemble them into a
[`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with an `otu_table` of per-sample abundances and the
representative sequence of each variant in the `refseq` slot. Three
denoising engines are available:

- `method = "vsearch"` (default):
  [vsearch](https://github.com/torognes/vsearch) UNOISE3 denoising
  (`--cluster_unoise`), producing zOTUs. Works on fasta, fastq and
  `dna_seq` inputs.

- `method = "swarm"`: [SWARM](https://github.com/torognes/swarm) with
  `d = 1` (fastidious by default). Works on fasta, fastq and `dna_seq`
  inputs.

- `method = "dada2"`: the classical
  [dada2](https://benjjneb.github.io/dada2/) denoising pipeline
  ([`dada2::dada()`](https://rdrr.io/pkg/dada2/man/dada.html) with
  `selfConsist = TRUE` by default, overridable through `...`), followed
  by de-novo chimera removal
  ([`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html)).
  Requires fastq files (quality scores drive the dada2 error model);
  reads containing non-ACGT bases are discarded with a message.

For fasta/fastq input, reads are first dereplicated per sample with
vsearch, then pooled and dereplicated again so that abundances per
sample can be rebuilt after denoising (see the details of
[`cluster_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cluster_reads.md)).
Use
[`cluster_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cluster_reads.md)
instead of this function to cluster reads into similarity-threshold OTUs
(vsearch `--cluster_size` or SWARM with `d > 1`).

## Usage

``` r
denoised_reads(
  path_to_fastx = NULL,
  dna_seq = NULL,
  sample_name = "sample_1",
  method = c("vsearch", "swarm", "dada2"),
  unoise_alpha = 2,
  d = 1,
  fastidious = TRUE,
  minsize = 1,
  nproc = 1,
  remove_chimeras = TRUE,
  vsearchpath = find_vsearch(),
  swarmpath = "swarm",
  swarm_args = "",
  keep_temporary_files = FALSE,
  ...
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
  abundance as a `;size=N` suffix in their names. Not available for
  `method = "dada2"` (no quality scores). Set either `path_to_fastx` or
  `dna_seq`, not both.

- sample_name:

  (default: `"sample_1"`) Sample name used when `dna_seq` is set.

- method:

  (default: `"vsearch"`) Denoising engine: `"vsearch"` (UNOISE3),
  `"swarm"` or `"dada2"`.

- unoise_alpha:

  (numeric, default 2) The `alpha` parameter of vsearch UNOISE
  (trade-off between sensitivity and specificity).

- d:

  (default: 1) SWARM resolution. `d = 1` is the denoising resolution;
  larger values cluster at a lower resolution and trigger a warning
  pointing to
  [`cluster_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cluster_reads.md).

- fastidious:

  (logical, default: TRUE) Use SWARM fastidious mode. Ignored when
  `d != 1`.

- minsize:

  (default: 1) Discard dereplicated sequences with an abundance strictly
  lower than `minsize` before denoising. Also passed to vsearch UNOISE
  (`--minsize`).

- nproc:

  (default: 1) Number of threads (SWARM and dada2 chimera removal).

- remove_chimeras:

  (logical, default: TRUE) For `method = "dada2"`, remove chimeras de
  novo with
  [`dada2::removeBimeraDenovo()`](https://rdrr.io/pkg/dada2/man/removeBimeraDenovo.html)
  after denoising.

- vsearchpath:

  (default:
  [`find_vsearch()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/find_vsearch.md))
  Path to the vsearch binary.

- swarmpath:

  (default: `"swarm"`) Path to the swarm binary.

- swarm_args:

  (default: `""`) Extra arguments passed to swarm.

- keep_temporary_files:

  (logical, default: FALSE) Keep the dereplicated fasta files, the
  pooled fasta, the cluster table and the representative-sequences fasta
  in a temporary directory.

- ...:

  Additional arguments passed to
  [`dada2::dada()`](https://rdrr.io/pkg/dada2/man/dada.html) when
  `method = "dada2"` (e.g. `err`, `selfConsist`, `OMEGA_A`, `pool`).

## Value

A
[`phyloseq::phyloseq-class()`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
object with one row (taxon) per ASV/zOTU named `ASV_1`, `ASV_2`, ...
ordered by decreasing total abundance, one column (sample) per input
file, the representative sequence of each variant in the `refseq` slot,
and (for file input) the source file paths in the `sample_data` slot.

## Details

vsearch is required for all methods (per-sample dereplication and, for
`method = "vsearch"`, UNOISE denoising); swarm is additionally required
for `method = "swarm"`. Use
[`is_vsearch_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_vsearch_installed.md)
and
[`is_swarm_installed()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/is_swarm_installed.md)
to check availability. This function is mainly a wrapper of the work of
others: please cite [vsearch](https://github.com/torognes/vsearch)
(<https://pubmed.ncbi.nlm.nih.gov/27781170>),
[SWARM](https://github.com/torognes/swarm)
([doi:10.1093/bioinformatics/btab493](https://doi.org/10.1093/bioinformatics/btab493)
) and/or [dada2](https://doi.org/10.1038/nmeth.3869) depending on the
method used.

## See also

[`cluster_reads()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cluster_reads.md),
[`asv2otu()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/postcluster_pq.md),
[`swarm_clustering()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/swarm_clustering.md),
[`vsearch_clustering()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/vsearch_clustering.md),
[`chimera_removal_vs()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/chimera_removal_vs.md)

## Author

Adrien Taudière

## Examples

``` r
f1 <- system.file("extdata", "ex_R1_001.fasta", package = "MiscMetabar")
f2 <- system.file("extdata", "ex_R1_002.fasta", package = "MiscMetabar")
(pq_zotu <- denoised_reads(path_to_fastx = c(f1, f2)))
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 50 taxa and 2 samples ]
#> sample_data() Sample Data:       [ 2 samples by 1 sample variables ]
#> refseq()      DNAStringSet:      [ 50 reference sequences ]

denoised_reads(path_to_fastx = c(f1, f2), method = "swarm")
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 50 taxa and 2 samples ]
#> sample_data() Sample Data:       [ 2 samples by 1 sample variables ]
#> refseq()      DNAStringSet:      [ 50 reference sequences ]
# \donttest{
if (requireNamespace("dada2")) {
  fq1 <- system.file("extdata", "ex_R1_001.fastq.gz", package = "MiscMetabar")
  fq2 <- system.file("extdata", "ex.fastq", package = "MiscMetabar")
  denoised_reads(path_to_fastx = c(fq1, fq2), method = "dada2")
}
#> "ex": discarding 1 unique sequence containing non-ACGT bases before dada2
#> denoising.
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 2 taxa and 2 samples ]
#> sample_data() Sample Data:       [ 2 samples by 1 sample variables ]
#> refseq()      DNAStringSet:      [ 2 reference sequences ]
# }
```
