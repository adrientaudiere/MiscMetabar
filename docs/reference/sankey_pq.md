# Sankey plot of [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html) object

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Graphical representation of distribution of taxa across Taxonomy and
(optionnaly a factor).

## Usage

``` r
sankey_pq(
  physeq = NULL,
  fact = NULL,
  taxa = 1:4,
  add_nb_seq = FALSE,
  min_prop_tax = 0,
  tax2remove = NULL,
  units = NULL,
  symbol2sub = c("\\.", "-"),
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  Name of the factor to cluster samples by modalities. Need to be in
  `physeq@sam_data`.

- taxa:

  a vector of taxonomic rank to plot

- add_nb_seq:

  Represent the number of sequences or the number of OTUs (add_nb_seq =
  FALSE). Note that plotting the number of sequences is slower.

- min_prop_tax:

  (default: 0) The minimum proportion for taxa to be plotted.
  EXPERIMENTAL. For the moment each links below the min.prop. tax is
  discard from the sankey network resulting in sometimes weird plot.

- tax2remove:

  a vector of taxonomic groups to remove from the analysis (e.g.
  `c('Incertae sedis', 'unidentified')`)

- units:

  character string describing physical units (if any) for Value

- symbol2sub:

  (default: c('\\', '-')) vector of symbol to delete in the taxonomy

- ...:

  Additional arguments passed on to
  [`sankeyNetwork`](https://rdrr.io/pkg/networkD3/man/sankeyNetwork.html)

## Value

A
[`sankeyNetwork`](https://rdrr.io/pkg/networkD3/man/sankeyNetwork.html)
plot representing the taxonomic distribution of OTUs or sequences. If
`fact` is set, represent the distribution of the last taxonomic level in
the modalities of `fact`

## See also

[`sankeyNetwork`](https://rdrr.io/pkg/networkD3/man/sankeyNetwork.html),
[`ggaluv_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/ggaluv_pq.md)

## Author

Adrien Taudière

## Examples

``` r
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
if (requireNamespace("networkD3")) {
  sankey_pq(GP, fact = "SampleType")
}
#> Loading required namespace: networkD3
#> Warning: NAs introduced by coercion

{"x":{"links":{"source":[0,0,1,2,2,2,1,1,2,1,1,3,4,5,6,7,8,8,9,4,10,3,11,12,13,14,15,11,12,13,16,14,15,11,12,13,19,14,15,11,12,13,16,14,15,17,11,12,20,13,19,14,21,11,12,13,18,19,14,15,11,12,13,14,15],"target":[1,2,3,6,7,8,9,4,5,10,22,17,11,12,20,13,18,16,19,14,21,15,23,23,23,23,23,24,24,24,24,24,24,26,26,26,26,26,26,27,27,27,27,27,27,29,29,29,29,29,29,29,29,30,30,30,30,30,30,30,31,31,31,31,31],"value":[106,102,27,7,27,24,8,57,44,3,8,7,37,44,7,27,12,12,5,20,1,18,5,4,8,4,2,5,6,6,1,3,11,8,7,3,1,11,4,21,19,7,1,2,8,1,4,3,1,1,1,12,1,7,9,3,2,4,19,1,6,9,2,3,4]},"nodes":{"name":["Archaea","Crenarchaeota","Euryarchaeota","C2","Thaumarchaeota","Thermoplasmata","Halobacteria","Methanobacteria","Methanomicrobia","SdNA","Thermoprotei","Cenarchaeales","E2","Methanobacteriales","Nitrososphaerales","pGrfC26","Methanosarcinales","B10","Methanomicrobiales","NRPJ","Halobacteriales","Sulfolobales","pMC2A209","FECES","FRESHWATER","FRESHWATER (CREEK)","MOCK","OCEAN","SEDIMENT (ESTUARY)","SKIN","SOIL","TONGUE"],"group":["Archaea","Crenarchaeota","Euryarchaeota","C2","Thaumarchaeota","Thermoplasmata","Halobacteria","Methanobacteria","Methanomicrobia","SdNA","Thermoprotei","Cenarchaeales","E2","Methanobacteriales","Nitrososphaerales","pGrfC26","Methanosarcinales","B10","Methanomicrobiales","NRPJ","Halobacteriales","Sulfolobales","pMC2A209","FECES","FRESHWATER","FRESHWATER (CREEK)","MOCK","OCEAN","SEDIMENT (ESTUARY)","SKIN","SOIL","TONGUE"]},"options":{"NodeID":"name","NodeGroup":"name","LinkGroup":null,"colourScale":"d3.scaleOrdinal(d3.schemeCategory20);","fontSize":7,"fontFamily":null,"nodeWidth":15,"nodePadding":10,"units":"OTUs","margin":{"top":null,"right":null,"bottom":null,"left":null},"iterations":32,"sinksRight":true}},"evals":[],"jsHooks":[]}# \donttest{
if (requireNamespace("networkD3")) {
  sankey_pq(GP, taxa = 1:4, min_prop_tax = 0.01)
  sankey_pq(GP, taxa = 1:4, min_prop_tax = 0.01, add_nb_seq = TRUE)
}

{"x":{"links":{"source":[0,0,1,2,1,2,3,4,5,3,6],"target":[1,2,6,5,3,4,7,8,9,10,11],"value":[79945,115653,22138,10377,56403,104475,25607,104475,10377,30796,21834]},"nodes":{"name":["Archaea","Crenarchaeota","Euryarchaeota","Thaumarchaeota","Thermoplasmata","Methanobacteria","C2","Cenarchaeales","E2","Methanobacteriales","Nitrososphaerales","pGrfC26"],"group":["Archaea","Crenarchaeota","Euryarchaeota","Thaumarchaeota","Thermoplasmata","Methanobacteria","C2","Cenarchaeales","E2","Methanobacteriales","Nitrososphaerales","pGrfC26"]},"options":{"NodeID":"name","NodeGroup":"name","LinkGroup":null,"colourScale":"d3.scaleOrdinal(d3.schemeCategory20);","fontSize":7,"fontFamily":null,"nodeWidth":15,"nodePadding":10,"units":"Sequences","margin":{"top":null,"right":null,"bottom":null,"left":null},"iterations":32,"sinksRight":true}},"evals":[],"jsHooks":[]}# }
```
