# Venn diagram of [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object using `ggVennDiagram::ggVennDiagram` function

[![lifecycle-stable](https://img.shields.io/badge/lifecycle-stable-green)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that you can use ggplot2 function to customize the plot for ex.
`+ scale_fill_distiller(palette = "BuPu", direction = 1)` and
`+ scale_x_continuous(expand = expansion(mult = 0.5))`. See examples.

## Usage

``` r
ggvenn_pq(
  physeq = NULL,
  fact = NULL,
  min_nb_seq = 0,
  taxonomic_rank = NULL,
  split_by = NULL,
  add_nb_samples = TRUE,
  add_nb_seq = FALSE,
  rarefy_before_merging = FALSE,
  rarefy_after_merging = FALSE,
  rngseed = FALSE,
  return_data_for_venn = FALSE,
  verbose = TRUE,
  type = "nb_taxa",
  na_remove = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor to cluster samples by modalities. Need
  to be in `physeq@sam_data`.

- min_nb_seq:

  minimum number of sequences by OTUs by samples to take into count this
  OTUs in this sample. For example, if min_nb_seq=2,each value of 2 or
  less in the OTU table will not count in the venn diagram

- taxonomic_rank:

  Name (or number) of a taxonomic rank to count. If set to Null (the
  default) the number of OTUs is counted.

- split_by:

  Split into multiple plot using variable split_by. The name of a
  variable must be present in `sam_data` slot of the physeq object.

- add_nb_samples:

  (logical, default TRUE) Add the number of samples to levels names

- add_nb_seq:

  (logical, default FALSE) Add the number of sequences to levels names

- rarefy_before_merging:

  Rarefy each sample before merging by the modalities of args `fact`.
  Use
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  function

- rarefy_after_merging:

  Rarefy each sample after merging by the modalities of args `fact`.

- rngseed:

  (Optional). A single integer value passed to
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html),
  which is used to fix a seed for reproducibly random number generation
  (in this case, reproducibly random subsampling). If set to FALSE, then
  no fiddling with the RNG seed is performed, and it is up to the user
  to appropriately call set.seed beforehand to achieve reproducible
  results. Default is FALSE.

- return_data_for_venn:

  (logical, default FALSE) If TRUE, the plot is not returned, but the
  resulting dataframe to plot with ggVennDiagram package is returned.

- verbose:

  (logical, default TRUE) If TRUE, prompt some messages.

- type:

  If "nb_taxa" (default), the number of taxa (ASV, OTU or taxonomic_rank
  if `taxonomic_rank` is not NULL) is used in plot. If "nb_seq", the
  number of sequences is plotted. `taxonomic_rank` is never used if type
  = "nb_seq".

- na_remove:

  (logical, default TRUE) If set to TRUE, remove samples with NA in the
  variables set in `fact` param

- ...:

  Other arguments for the
  [`ggVennDiagram::ggVennDiagram`](https://gaospecial.github.io/ggVennDiagram/reference/ggVennDiagram.html)
  function for ex. `category.names`.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot
representing Venn diagram of modalities of the argument `factor` or if
split_by is set a list of plots.

## See also

[`upset_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/upset_pq.md)

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("ggVennDiagram")) {
  ggvenn_pq(data_fungi_mini, fact = "Height")
}
#> 47 were discarded due to NA in variable fact
#> Taxa are now in columns.
#> Cleaning suppress 11 taxa and 0 samples.
#> Cleaning suppress 11 taxa and 0 samples.
#> Cleaning suppress 11 taxa and 0 samples.
```
