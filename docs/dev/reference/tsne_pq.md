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

## Details

t-SNE is a **local** dimensionality-reduction technique: it
preferentially preserves the neighbourhood structure of the original
samples. It is therefore well suited to *neighbourhood identification*,
*outlier identification* and *cluster identification*. However, the
distances between points, the distances between clusters, the density of
a cluster and the class separability read on the embedding are **not**
faithful to the original space and must not be interpreted as such (Jeon
et al., 2026). For these *global* questions, use a global technique
instead, such as PCA
([`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html)) or MDS/PCoA
([`phyloseq::ordinate()`](https://rdrr.io/pkg/phyloseq/man/ordinate.html),
[`plot_ordination_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/plot_ordination_pq.md)).
Recent techniques such as UMATO (Jeon et al., 2025) and PaCMAP (Wang et
al., 2021) seek a better compromise between local and global structure
preservation.

This function is mainly a wrapper of the work of others. Please make a
reference to
[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html) if you use
this function.

## References

Jeon, H., Park, J., Shin, S., & Seo, J. (2026). Stop Misusing t-SNE and
UMAP for Visual Analytics. *IEEE Transactions on Visualization and
Computer Graphics*.
[doi:10.48550/arXiv.2506.08725](https://doi.org/10.48550/arXiv.2506.08725)

van der Maaten, L., & Hinton, G. (2008). Visualizing Data using t-SNE.
*Journal of Machine Learning Research*, 9, 2579-2605.

## See also

[`Rtsne::Rtsne()`](https://rdrr.io/pkg/Rtsne/man/Rtsne.html),
[`umap_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/umap_pq.md),
[`plot_ordination_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/plot_ordination_pq.md),
[`phyloseq::plot_ordination()`](https://rdrr.io/pkg/phyloseq/man/plot_ordination.html)

## Author

Adrien Taudière

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
