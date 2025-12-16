# Build phylogenetic trees from refseq slot of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This function build tree phylogenetic tree and if nb_bootstrap is set,
it build also the 3 corresponding bootstrapped tree.

Default parameters are based on
[doi:10.12688/f1000research.8986.2](https://doi.org/10.12688/f1000research.8986.2)
and phangorn vignette [Estimating phylogenetic trees with
phangorn](https://klausvigo.github.io/phangorn/articles/Trees.html). You
should understand your data, especially the markers, before using this
function.

Note that phylogenetic reconstruction with markers used for
metabarcoding are not robust. You must verify the robustness of your
phylogenetic tree using taxonomic classification (see vignette [Tree
visualization](https://adrientaudiere.github.io/MiscMetabar/articles/tree_visualization.html))
and bootstrap or multi-tree visualization

## Usage

``` r
build_phytree_pq(
  physeq,
  nb_bootstrap = 0,
  model = "GTR",
  optInv = TRUE,
  optGamma = TRUE,
  rearrangement = "NNI",
  control = phangorn::pml.control(trace = 0),
  optNni = TRUE,
  multicore = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- nb_bootstrap:

  (default 0): If a positive number is set, the function also build 3
  bootstrapped trees using `nb_bootstrap` bootstrap samples

- model:

  allows to choose an amino acid models or nucleotide model, see
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  for more details

- optInv:

  Logical value indicating whether topology gets optimized (NNI). See
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  for more details

- optGamma:

  Logical value indicating whether gamma rate parameter gets optimized.
  See
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  for more details

- rearrangement:

  type of tree tree rearrangements to perform, one of "NNI",
  "stochastic" or "ratchet" see
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  for more details

- control:

  A list of parameters for controlling the fitting process. see
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  for more details

- optNni:

  Logical value indicating whether topology gets optimized (NNI). see
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  for more details

- multicore:

  (logical) whether models should estimated in parallel. see
  [`phangorn::bootstrap.pml()`](https://klausvigo.github.io/phangorn/reference/bootstrap.pml.html)
  for more details

- ...:

  Other params for be passed on to
  [`phangorn::optim.pml()`](https://klausvigo.github.io/phangorn/reference/pml.html)
  function

## Value

A list of phylogenetic tree

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to `phangorn` package if you use this function.

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("phangorn")) {
  set.seed(22)
  df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
  df_tree <- build_phytree_pq(df, nb_bootstrap = 2)
  plot(df_tree$UPGMA)
  phangorn::plotBS(df_tree$UPGMA, df_tree$UPGMA_bs, main = "UPGMA")
  plot(df_tree$NJ, "unrooted")
  plot(df_tree$ML)

  phangorn::plotBS(df_tree$ML$tree, df_tree$ML_bs, p = 20, frame = "circle")
  phangorn::plotBS(
    df_tree$ML$tree,
    df_tree$ML_bs,
    p = 20,
    frame = "circle",
    method = "TBE"
  )
  plot(phangorn::consensusNet(df_tree$ML_bs))
  plot(phangorn::consensusNet(df_tree$NJ_bs))
  ps_tree <- merge_phyloseq(df, df_tree$ML$tree)
}
#> Cleaning suppress 0 taxa (  ) and 6 sample(s) ( AD26-005-H_S10_MERGED.fastq.gz / CB8-019-H_S70_MERGED.fastq.gz / DY5-004-H_S97_MERGED.fastq.gz / N23-002-B_S130_MERGED.fastq.gz / NVABM0244-M_S137_MERGED.fastq.gz / T28-ABM602-B_S162_MERGED.fastq.gz ).
#> Number of non-matching ASV 0
#> Number of matching ASV 45
#> Number of filtered-out ASV 23
#> Number of kept ASV 22
#> Number of kept samples 131
#> Determining distance matrix based on shared 8-mers:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Clustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.01 secs
#> 
#> Aligning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.16 secs
#> 
#> Iteration 1 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0.01 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.11 secs
#> 
#> Iteration 2 of 2:
#> 
#> Determining distance matrix based on alignment:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Reclustering into groups by similarity:
#> ================================================================================
#> 
#> Time difference of 0 secs
#> 
#> Realigning Sequences:
#> ================================================================================
#> 
#> Time difference of 0.08 secs
#> 
#> Refining the alignment:
#> ================================================================================
#> 
#> Time difference of 0.03 secs
#> 
#> optimize edge weights:  -4496.258 --> -4357.83 
#> optimize edge weights:  -4357.83 --> -4357.827 
#> optimize topology:  -4357.827 --> -4357.827  NNI moves:  0 
#> optimize edge weights:  -4357.827 --> -4357.827 
#> optimize edge weights:  -4438.165 --> -4308.211 
#> optimize edge weights:  -4308.211 --> -4308.21 
#> optimize topology:  -4308.21 --> -4308.21  NNI moves:  0 
#> optimize edge weights:  -4308.21 --> -4308.21 








#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
# }
```
