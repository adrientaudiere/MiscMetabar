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
  neighbors (Default 15). See
  [`umap::umap.defaults()`](https://rdrr.io/pkg/umap/man/umap.defaults.html)
  or
  [`uwot::umap2()`](https://jlmelville.github.io/uwot/reference/umap2.html)
  for the list of parameters and default values.

## Value

A dataframe with samples informations and the x_umap and y_umap position

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to [`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html) if
you use this function.

## See also

[`umap::umap()`](https://rdrr.io/pkg/umap/man/umap.html),
[`tsne_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/tsne_pq.md),
[`phyloseq::plot_ordination()`](https://rdrr.io/pkg/phyloseq/man/plot_ordination.html)

## Author

Adrien Taudière

## Examples

``` r
library("umap")
df_umap <- umap_pq(data_fungi_mini)
#> Taxa are now in columns.
#> Taxa are now in rows.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Joining with `by = join_by(Sample)`
#> Joining with `by = join_by(Sample)`
ggplot(df_umap, aes(x = x_umap, y = y_umap, col = Height)) +
  geom_point(size = 2)


# \donttest{
library(patchwork)
physeq <- data_fungi_mini
res_tsne <- tsne_pq(data_fungi_mini)
df_umap_tsne <- df_umap
df_umap_tsne$x_tsne <- res_tsne$Y[, 1]
df_umap_tsne$y_tsne <- res_tsne$Y[, 2]
((ggplot(df_umap, aes(x = x_umap, y = y_umap, col = Height)) +
  geom_point(size = 2) +
  ggtitle("UMAP")) + (plot_ordination(physeq,
  ordination =
    ordinate(physeq, method = "PCoA", distance = "bray"), color = "Height"
)) +
  ggtitle("PCoA")) /
  ((ggplot(df_umap_tsne, aes(x = x_tsne, y = y_tsne, col = Height)) +
    geom_point(size = 2) +
    ggtitle("tsne")) +
    (plot_ordination(physeq,
      ordination = ordinate(physeq, method = "NMDS", distance = "bray"),
      color = "Height"
    ) +
      ggtitle("NMDS"))) +
  patchwork::plot_layout(guides = "collect")
#> Square root transformation
#> Wisconsin double standardization
#> Run 0 stress 0.1842417 
#> Run 1 stress 0.1829537 
#> ... New best solution
#> ... Procrustes: rmse 0.06810274  max resid 0.2633582 
#> Run 2 stress 0.1825839 
#> ... New best solution
#> ... Procrustes: rmse 0.07838589  max resid 0.3894694 
#> Run 3 stress 0.1875225 
#> Run 4 stress 0.186983 
#> Run 5 stress 0.1814762 
#> ... New best solution
#> ... Procrustes: rmse 0.05418813  max resid 0.2315781 
#> Run 6 stress 0.1855187 
#> Run 7 stress 0.1834722 
#> Run 8 stress 0.1826215 
#> Run 9 stress 0.1851733 
#> Run 10 stress 0.1866514 
#> Run 11 stress 0.1807717 
#> ... New best solution
#> ... Procrustes: rmse 0.06219345  max resid 0.3913131 
#> Run 12 stress 0.1861972 
#> Run 13 stress 0.1847311 
#> Run 14 stress 0.1839897 
#> Run 15 stress 0.1866369 
#> Run 16 stress 0.1856367 
#> Run 17 stress 0.1817974 
#> Run 18 stress 0.184205 
#> Run 19 stress 0.1861744 
#> Run 20 stress 0.1832293 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     17: no. of iterations >= maxit
#>      3: stress ratio > sratmax


df_uwot <- umap_pq(data_fungi_mini, pkg = "uwot")
#> Taxa are now in columns.
#> Taxa are now in rows.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Joining with `by = join_by(Sample)`
#> Joining with `by = join_by(Sample)`

(ggplot(df_umap, aes(x = x_umap, y = y_umap, col = Height)) +
  geom_point(size = 2) +
  ggtitle("umap::umap")) /
  (ggplot(df_uwot, aes(x = x_umap, y = y_umap, col = Height)) +
    geom_point(size = 2) +
    ggtitle("uwot::umap2"))

# }
```
