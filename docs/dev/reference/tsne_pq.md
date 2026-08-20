# Compute tSNE position of samples from a phyloseq object

Compute tSNE position of samples from a phyloseq object

## Usage

``` r
tsne_pq(physeq, method = "bray", dims = 2, theta = 0, perplexity = 30, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- method:

  A method to calculate distance using
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function

- dims:

  (Int) Output dimensionality (default: 2)

- theta:

  (Numeric) Speed/accuracy trade-off (increase for less accuracy), set
  to 0.0 for exact TSNE (default: 0.0 see details in the man page of
  [`Rtsne::Rtsne`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html)).

- perplexity:

  (Numeric) Perplexity parameter (should not be bigger than 3 \*
  perplexity \< nrow(X) - 1, see details in the man page of
  [`Rtsne::Rtsne`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html))

- ...:

  Additional arguments passed on to
  [`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html)

## Value

A dataframe with samples informations and the x_tsne and y_tsne position
(or `tsne_1`, `tsne_2`, ... columns when `dims != 2`).

## See also

[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html),
[`umap_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/umap_pq.md),
[`phyloseq::plot_ordination()`](https://rdrr.io/pkg/phyloseq/man/plot_ordination.html)

## Author

Adrien Taudiere

## Examples

``` r
if (requireNamespace("Rtsne")) {
  df_tsne <- tsne_pq(data_fungi_mini)
  ggplot(df_tsne, aes(x = x_tsne, y = y_tsne, col = Height)) +
    geom_point(size = 2)
}
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Joining with `by = join_by(Sample)`
#> Joining with `by = join_by(Sample)`
```
