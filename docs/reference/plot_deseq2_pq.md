# Plot DESeq2 results for a phyloseq or a DESeq2 object.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of DESeq2 analysis.

## Usage

``` r
plot_deseq2_pq(
  data,
  contrast = NULL,
  tax_table = NULL,
  pval = 0.05,
  taxolev = "Genus",
  select_taxa = NULL,
  color_tax = "Phylum",
  tax_depth = NULL,
  verbose = TRUE,
  jitter_width = 0.1,
  ...
)
```

## Arguments

- data:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  or a
  [`DESeqDataSet-class`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html)
  object.

- contrast:

  (required) contrast specifies what comparison to extract from the
  object to build a results table. See
  [`results`](https://rdrr.io/pkg/DESeq2/man/results.html) man page for
  more details.

- tax_table:

  Required if data is a
  [`DESeqDataSet-class`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html)
  object. The taxonomic table used to find the `taxa` and `color_taxa`
  arguments. If data is a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object, data@tax_table is used.

- pval:

  (default: 0.05) the significance cutoff used for optimizing the
  independent filtering. If the adjusted p-value cutoff (FDR) will be a
  value other than 0.05, pval should be set to that value.

- taxolev:

  taxonomic level of interest

- select_taxa:

  Either the name of the taxa (in the form of
  [`DESeq2::results()`](https://rdrr.io/pkg/DESeq2/man/results.html)) or
  a logical vector (length of the results from
  [`DESeq2::results()`](https://rdrr.io/pkg/DESeq2/man/results.html)) to
  select taxa to plot.

- color_tax:

  taxonomic level used for color or a color vector.

- tax_depth:

  Taxonomic depth to test for differential distribution among contrast.
  If Null the analysis is done at the OTU (i.e. Species) level. If not
  Null, data need to be a column name in the `tax_table` slot of the
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object.

- verbose:

  whether the function print some information during the computation

- jitter_width:

  width for the jitter positioning

- ...:

  Additional arguments passed on to
  [`DESeq`](https://rdrr.io/pkg/DESeq2/man/DESeq.html) or
  [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot
representing DESeq2 results

## Details

Please cite `DESeq2` package if you use chis function.

## See also

[`DESeq`](https://rdrr.io/pkg/DESeq2/man/DESeq.html)

[`results`](https://rdrr.io/pkg/DESeq2/man/results.html)

[`plot_edgeR_pq`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_edgeR_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{

data("GlobalPatterns", package = "phyloseq")
GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
GP <- subset_samples(GP, SampleType %in% c("Soil", "Skin"))
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
if (requireNamespace("DESeq2")) {
  res <- DESeq2::DESeq(phyloseq_to_deseq2(GP, ~SampleType),
    test = "Wald", fitType = "local"
  )
  plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"),
    tax_table = GP@tax_table, color_tax = "Kingdom"
  )
  plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"),
    tax_table = GP@tax_table, color_tax = "Kingdom",
    pval = 0.7
  )
  plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"),
    tax_table = GP@tax_table, color_tax = "Class",
    select_taxa = c("522457", "271582")
  )
}
#> Loading required namespace: DESeq2
#> converting counts to integer mode
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
#> None taxa present significant distribution pattern through
#>               contrast.
#> [1] "None taxa present significant distribution pattern through\n             contrast."
# }
```
