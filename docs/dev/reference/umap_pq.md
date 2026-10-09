# Computes a manifold approximation and projection (UMAP) for phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

https://journals.asm.org/doi/full/10.1128/msystems.00691-21

## Usage

``` r
umap_pq(physeq, pkg = "umap", ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- pkg:

  Which R packages to use, either "umap" or "uwot".

- ...:

  Additional arguments passed on to
  [`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html) or
  [`uwot::umap2()`](https://jlmelville.github.io/uwot/reference/umap2.html)
  function. For example `n_neighbors` set the number of nearest
  neighbors (Default 15). `n_neighbors` is capped to
  `nsamples(physeq) - 1` (with a minimum of 2) so that small phyloseq
  objects run, as in
  [`ggplotpq::dr_plot_pq()`](https://adrientaudiere.github.io/ggplotpq/reference/dr_plot_pq.html).
  See
  [`umap::umap.defaults()`](https://rdrr.io/pkg/umap/man/umap.defaults.html)
  or
  [`uwot::umap2()`](https://jlmelville.github.io/uwot/reference/umap2.html)
  for the list of parameters and default values.

## Value

A dataframe with samples informations and the x_umap and y_umap position

## Details

UMAP is a **local** dimensionality-reduction technique: it
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
reference to [`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html) if
you use this function.

## References

Jeon, H., Park, J., Shin, S., & Seo, J. (2026). Stop Misusing t-SNE and
UMAP for Visual Analytics. *IEEE Transactions on Visualization and
Computer Graphics*.
[doi:10.48550/arXiv.2506.08725](https://doi.org/10.48550/arXiv.2506.08725)

Jeon, H., Ko, K., Lee, S., Hyun, J., Yang, T., Go, G., Jo, J., & Seo, J.
(2025). UMATO: Bridging Local and Global Structures for Reliable Visual
Analytics with Dimensionality Reduction. *IEEE Transactions on
Visualization and Computer Graphics*.
[doi:10.1109/TVCG.2025.3602735](https://doi.org/10.1109/TVCG.2025.3602735)

McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold
Approximation and Projection for Dimension Reduction.
[doi:10.48550/arXiv.1802.03426](https://doi.org/10.48550/arXiv.1802.03426)

Wang, Y., Huang, H., Rudin, C., & Shaposhnik, Y. (2021). Understanding
How Dimension Reduction Tools Work: An Empirical Approach to Deciphering
t-SNE, UMAP, TriMap, and PaCMAP for Data Visualization. *Journal of
Machine Learning Research*, 22(201), 1-73.
[doi:10.48550/arXiv.2012.04456](https://doi.org/10.48550/arXiv.2012.04456)

## See also

[`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html),
[`tsne_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/tsne_pq.md),
[`phyloseq::plot_ordination()`](https://rdrr.io/pkg/phyloseq/man/plot_ordination.html)

## Author

Adrien Taudière

## Examples

``` r
library("umap")
data_f <- prune_samples(
  sample_names(data_fungi_mini)[1:20],
  data_fungi_mini
)
df_umap <- umap_pq(data_f, n_neighbors = 3)
#> Taxa are now in columns.
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Joining with `by = join_by(Sample)`
ggplot(df_umap, aes(x = x_umap, y = y_umap, col = Height)) +
  geom_point(size = 2)


if (FALSE) { # \dontrun{
df_uwot <- umap_pq(data_fungi_mini, pkg = "uwot")
library(patchwork)
physeq <- data_fungi_mini
df_umap <- umap_pq(physeq, n_neighbors = 3)
df_umap_tsne <- tsne_pq(data_fungi_mini)
((ggplot(df_umap, aes(x = x_umap, y = y_umap, col = Height)) +
  geom_point(size = 2) +
  ggtitle("UMAP")) +
  (plot_ordination(physeq,
    ordination = ordinate(physeq, method = "PCoA", distance = "bray"),
    color = "Height"
  ) + ggtitle("PCoA"))) /
  ((ggplot(df_umap_tsne, aes(x = x_tsne, y = y_tsne, col = Height)) +
    geom_point(size = 2) +
    ggtitle("tsne")) +
    (plot_ordination(physeq,
      ordination = ordinate(physeq, method = "NMDS", distance = "bray"),
      color = "Height"
    ) + ggtitle("NMDS"))) +
  patchwork::plot_layout(guides = "collect")

(ggplot(df_umap, aes(x = x_umap, y = y_umap, col = Height)) +
  geom_point(size = 2) +
  ggtitle("umap::umap")) /
  (ggplot(df_uwot, aes(x = x_umap, y = y_umap, col = Height)) +
    geom_point(size = 2) +
    ggtitle("uwot::umap2"))
} # }
```
