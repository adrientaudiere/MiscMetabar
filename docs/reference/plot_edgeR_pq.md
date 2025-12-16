# Plot edgeR results for a phyloseq or a edgeR object.

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of edgeR result.

## Usage

``` r
plot_edgeR_pq(
  physeq,
  contrast = NULL,
  pval = 0.05,
  taxolev = "Genus",
  color_tax = "Phylum",
  verbose = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- contrast:

  (required) This argument specifies what comparison to extract from the
  object to build a results table. See
  [`results`](https://rdrr.io/pkg/DESeq2/man/results.html) man page for
  more details.

- pval:

  (default: 0.05) the significance cutoff used for optimizing the
  independent filtering. If the adjusted p-value cutoff (FDR) will be a
  value other than 0.05, pval should be set to that value.

- taxolev:

  taxonomic level of interest

- color_tax:

  taxonomic level used for color assignation

- verbose:

  (logical): whether the function print some information during the
  computation

- ...:

  Additional arguments passed on to
  [`exactTest`](https://rdrr.io/pkg/edgeR/man/exactTest.html) or
  [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot
representing edgeR results

## See also

[`exactTest`](https://rdrr.io/pkg/edgeR/man/exactTest.html)

[`plot_deseq2_pq`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_deseq2_pq.md)

## Author

Adrien Taudière

## Examples

``` r
data("GlobalPatterns", package = "phyloseq")
GP_archae <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
# \donttest{
if (requireNamespace("edgeR")) {
  plot_edgeR_pq(GP_archae, c("SampleType", "Soil", "Feces"),
    color_tax = "Kingdom"
  )

  plot_edgeR_pq(GP_archae, c("SampleType", "Soil", "Feces"),
    taxolev = "Class", color_tax = "Kingdom"
  )
}
#> Conversion to edgeR format
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Perform edgeR binary test
#> Conversion to edgeR format
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Perform edgeR binary test

# }
```
