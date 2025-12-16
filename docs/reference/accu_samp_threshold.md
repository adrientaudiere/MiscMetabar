# Compute the number of sequence to obtain a given proportion of ASV in accumulation curves

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Note that as most bioinformatic pipeline discard singleton, accumulation
curves from metabarcoding cannot be interpreted in the same way as with
conventional biodiversity sampling techniques.

## Usage

``` r
accu_samp_threshold(res_accuplot, threshold = 0.95)
```

## Arguments

- res_accuplot:

  the result of the function accu_plot()

- threshold:

  the proportion of ASV to obtain in each samples

## Value

a value for each sample of the number of sequences needed to obtain
`threshold` proportion of the ASV

## See also

[`accu_plot()`](https://adrientaudiere.github.io/MiscMetabar/reference/accu_plot.md)

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
GP <- rarefy_even_depth(subset_samples_pq(GP, sample_sums(GP) > 3000))
#> You set `rngseed` to FALSE. Make sure you've set & recorded
#>  the random seed of your session for reproducibility.
#> See `?set.seed`
#> ...
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
#> 58OTUs were removed because they are no longer 
#> present in any sample after random subsampling
#> ...
#> Found more than one class "phylo" in cache; using the first, from namespace 'phyloseq'
#> Also defined by ‘RNeXML’
p <- accu_plot(GP, "SampleType", add_nb_seq = TRUE, by.fact = TRUE, step = 10)

val_threshold <- accu_samp_threshold(p)

summary(val_threshold)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>    2891    6154    7361    7006    8214   10411 

##'  Plot the number of sequences needed to accumulate 0.95% of ASV in 50%, 75%
##'  and 100% of samples
p + geom_vline(xintercept = quantile(val_threshold, probs = c(0.50, 0.75, 1)))
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_ribbon()`).
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_line()`).

# }
```
