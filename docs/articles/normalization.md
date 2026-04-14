# Transformation and normalisation of phyloseq objects

``` r
library(MiscMetabar)
data(data_fungi_mini)
```

Choosing how to normalise or transform a count table is one of the most
consequential (and discussed) decisions in metabarcoding analysis.
Unequal sequencing depths across samples introduce a systematic bias
that can distort every downstream step — alpha diversity, ordination,
differential abundance testing. Yet every normalisation method makes its
own assumptions and trade-offs (McMurdie and Holmes 2014).

This article surveys the methods available in **MiscMetabar**, shows
their effect on `data_fungi_mini`, and offers a decision guide to help
you pick the right approach.

------------------------------------------------------------------------

## Why does it matter?

Amplicon sequencing returns *compositional* data: only the relative
proportions of taxa are observed, not their absolute abundances (Gloor
et al. 2017; Quinn et al. 2018). On top of that, libraries differ in
total read count by one to two orders of magnitude, even within the same
experiment. Applying diversity or ordination methods to raw counts
therefore conflates biological signal with technical depth variation.

``` r
ggplot(data=tibble(x=sample_sums(data_fungi_mini)),aes(x=x)) + geom_histogram(color="black") + scale_x_log10() + labs(x="Number of sequences per samples")
```

![Large variation in sequencing depth across samples in
'data_fungi_mini'.](normalization_files/figure-html/depth-variation-1.png)

Large variation in sequencing depth across samples in ‘data_fungi_mini’.

------------------------------------------------------------------------

## Overview of methods

| Family | Function | Key property |
|----|----|----|
| Presence/absence | [`as_binary_otu_table()`](https://adrientaudiere.github.io/MiscMetabar/reference/as_binary_otu_table.md) | Ignores abundance entirely |
| Proportions (TSS) | `transform_pq(method="tss")` | Divides by library size |
| TSS + log | [`normalize_prop_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/normalize_prop_pq.md) | TSS × constant, then log |
| Hellinger | `transform_pq(method="hellinger")` | Square-root of proportions |
| Centred log-ratio | `transform_pq(method="clr")` | Compositionally coherent |
| Robust CLR | `transform_pq(method="rclr")` | CLR robust to zeros |
| Log(1+x) | `transform_pq(method="log1p")` | Simple variance stabilisation |
| Z-score | `transform_pq(method="z")` | Per-taxon standardisation |
| Rank | `transform_pq(method="rank")` | Non-parametric, outlier robust |
| Rarefaction | [`rarefy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rarefy_pq.md) | Subsampling to equal depth |
| SRS | [`srs_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/srs_pq.md) | Rank-preserving subsampling (Heidrich et al. 2021) |
| GMPR | [`gmpr_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/gmpr_pq.md) | Pairwise ratio geometric mean (Chen et al. 2018) |
| CSS | [`css_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/css_pq.md) | Cumulative sum scaling (Paulson et al. 2013) |
| TMM | [`tmm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/tmm_pq.md) | Trimmed mean of M-values (Robinson and Oshlack 2010) |
| VST | [`vst_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/vst_pq.md) | Variance-stabilising (DESeq2) (Love, Huber, and Anders 2014) |
| Depth residuals | [`mcknight_residuals_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/mcknight_residuals_pq.md) | Log-log regression residuals (McKnight et al. 2019) |

------------------------------------------------------------------------

## Exploring the transformations

### `transform_pq()` — a unified interface via `vegan::decostand`

[`transform_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/transform_pq.md)
wraps
[`vegan::decostand()`](https://vegandevs.github.io/vegan/reference/decostand.html)
and handles the `taxa_are_rows` orientation automatically.

``` r
data_tss <- transform_pq(data_fungi_mini, method = "tss")
data_hell <- transform_pq(data_fungi_mini, method = "hellinger")
data_clr <- transform_pq(data_fungi_mini, method = "clr")  # pseudocount = 1 by default
data_log1p <- transform_pq(data_fungi_mini, method = "log1p")

# Sample sums after TSS: all 1
round(range(sample_sums(data_tss)), 4)
#> [1] 1 1
```

**Hellinger** and **CLR** are particularly recommended for ordination
(Legendre and Gallagher 2001; Gloor et al. 2017). Hellinger distance on
the transformed table equals the Hellinger distance without
transformation, which behaves better than Bray-Curtis on sparse data.
CLR makes the data compositionally coherent but requires replacing zeros
(vegan uses a pseudocount internally).

### `normalize_prop_pq()` — TSS with log scaling (McKnight 2019)

Multiplies relative abundances by a constant (default 10 000) then
applies a log transformation. This is the “rarefaction-free”
normalisation proposed by McKnight et al. (2019).

``` r
data_norm <- normalize_prop_pq(data_fungi_mini)
# All sample sums are equal after log transform
round(range(sample_sums(data_norm)), 2)
#> [1] 13.29 78.11
```

### `rarefy_pq()` — rarefaction with optional averaging

Rarefaction subsamples all libraries to the same depth, discarding
reads. It remains the most widely used method despite debates about data
loss (McMurdie and Holmes 2014; McKnight et al. 2019).
[`rarefy_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rarefy_pq.md)
allows averaging over `n` random subsamplings to reduce stochasticity.

``` r
data_rar1 <- rarefy_pq(data_fungi_mini, seed = 1)
data_rar10 <- rarefy_pq(data_fungi_mini, n = 10, seed = 1)

# Single rarefaction: all depths equal
unique(sample_sums(data_rar1))
#> [1] 1
# Averaged: sample sums close but not integer
round(range(sample_sums(data_rar10)), 1)
#> [1] 1 1
```

### `srs_pq()` — Scaling with Ranked Subsampling

SRS (Heidrich et al. 2021) scales to a common library size while
**preserving the rank order** of OTU abundances, avoiding the randomness
of rarefaction and the data loss associated with it.

``` r
data_srs <- srs_pq(data_fungi_mini)
round(range(sample_sums(data_srs)), 0)
#> [1] 1 1
```

### `gmpr_pq()` — Geometric Mean of Pairwise Ratios

GMPR (Chen et al. 2018) estimates size factors from the geometric mean
of median pairwise ratios between samples, making it robust to the high
proportion of zeros in microbial OTU tables.

``` r
data_gmpr <- gmpr_pq(data_fungi_mini)
round(range(sample_sums(data_gmpr)), 0)
#> [1]     1 61091
```

### `css_pq()` — Cumulative Sum Scaling

CSS (Paulson et al. 2013) normalises each sample by the cumulative sum
of counts up to a data-driven percentile, rather than the total. It is
less sensitive to a few very abundant OTUs.

``` r
data_css <- css_pq(data_fungi_mini)
```

### `tmm_pq()` — Trimmed Mean of M-values

TMM (Robinson and Oshlack 2010), borrowed from RNA-seq, computes
library-specific normalisation factors by trimming the distribution of
log-fold changes between each sample and a reference.

``` r
data_tmm <- tmm_pq(data_fungi_mini)
```

### `vst_pq()` — Variance Stabilising Transformation

VST (Love, Huber, and Anders 2014) fits a negative-binomial model to
stabilise the mean-variance relationship, making counts more suitable
for methods that assume homoscedastic data.

``` r
data_vst <- vst_pq(data_fungi_mini)
```

### `mcknight_residuals_pq()` — depth-robust alpha diversity

Instead of normalising the count table, this function computes residuals
from a regression of log-richness on log-depth and stores them in
`sample_data`. The residuals represent richness corrected for depth
variation — a depth-robust alpha diversity metric proposed by McKnight
et al. (2019) and used in large-scale soil fungal surveys (Mikryukov et
al. 2023).

``` r
data_res <- mcknight_residuals_pq(data_fungi_mini)
head(sample_data(data_res)$mcknight_residuals)
#> [1]  0.3273090 -0.1430749  0.5681529  0.7310108 -0.7698491  0.4374109
```

------------------------------------------------------------------------

## Comparing methods via ordination

A PCoA on Bray-Curtis dissimilarity illustrates how strongly the choice
of normalisation shapes the ordination landscape.

``` r
ord_raw <- phyloseq::ordinate(data_fungi_mini, method = "PCoA", distance = "bray")
ord_tss <- phyloseq::ordinate(data_tss, method = "PCoA", distance = "bray")
ord_hell <- phyloseq::ordinate(data_hell, method = "PCoA", distance = "bray")
ord_log1p <- phyloseq::ordinate(data_log1p, method = "PCoA", distance = "bray")

p_raw <- phyloseq::plot_ordination(
  data_fungi_mini, ord_raw, color = "Height"
) + ggplot2::ggtitle("Raw counts")

p_tss <- phyloseq::plot_ordination(
  data_tss, ord_tss, color = "Height"
) + ggplot2::ggtitle("TSS")

p_hell <- phyloseq::plot_ordination(
  data_hell, ord_hell, color = "Height"
) + ggplot2::ggtitle("Hellinger")

p_log1p <- phyloseq::plot_ordination(
  data_log1p, ord_log1p, color = "Height"
) + ggplot2::ggtitle("log1p")

patchwork::wrap_plots(p_raw, p_tss, p_hell, p_log1p, ncol = 2) +
  patchwork::plot_layout(guides = "collect")
```

![PCoA (Bray-Curtis) on raw counts, TSS, Hellinger, and
log1p-transformed data. Colour encodes the 'Height'
variable.](normalization_files/figure-html/ordination-comparison-1.png)

PCoA (Bray-Curtis) on raw counts, TSS, Hellinger, and log1p-transformed
data. Colour encodes the ‘Height’ variable.

------------------------------------------------------------------------

## Model-level depth correction (no table transformation)

Some MiscMetabar functions correct for sequencing depth **without**
transforming the OTU table, by incorporating library size directly into
the statistical model. These are complementary to — not replacements for
— the table-level normalisations above. If you transform the table, you
should not use these model-level corrections (set
correction_for_sample_size = FALSE).

### `adonis_pq()` — PERMANOVA with library-size covariate

When `correction_for_sample_size = TRUE`,
[`adonis_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/adonis_pq.md)
prepends `sample_size` (i.e. `sample_sums(physeq)`) to the PERMANOVA
formula, as recommended by Weiss et al. (2017):

    sample_size + Biological_Effect ~ distance_matrix

Library-size variation is partitioned out as the first term, so the
remaining variance attributed to the biological factor is not confounded
by depth. The raw count table is used as-is; only the formula changes.

### `hill_pq()` / `hill_tuckey_pq()` — alpha diversity with sqrt-depth residuals

When `correction_for_sample_size = TRUE` (the default), these functions:

1.  Compute Hill numbers from the **raw** OTU table.
2.  Fit `lm(hill_values ~ sqrt(read_numbers))` per diversity order.
3.  Use the **residuals** as the response in Tukey ANOVA.

The square-root of read count linearises the typical depth–richness
relationship; the residuals represent diversity corrected for depth
without discarding any data. Note that the Hill-number values
**displayed** are always computed from raw counts — only the
*statistical test* uses residuals.

This is conceptually similar to
[`mcknight_residuals_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/mcknight_residuals_pq.md)
but integrated into the plotting workflow rather than stored as a new
column in `sam_data`.

------------------------------------------------------------------------

## References

Chen, Li, James Reeve, Lunchao Zhang, Shengbing Huang, Xuefeng Wang, and
Jun Chen. 2018. “GMPR: A Robust Normalization Method for Zero-Inflated
Count Data with Application to Microbiome Sequencing Data.” *PeerJ* 6:
e4600. <https://doi.org/10.7717/peerj.4600>.

Gloor, Gregory B, Jean M Macklaim, Vera Pawlowsky-Glahn, and Juan J
Egozcue. 2017. “Microbiome Datasets Are Compositional: And This Is Not
Optional.” *Frontiers in Microbiology* 8: 2224.
<https://doi.org/10.3389/fmicb.2017.02224>.

Heidrich, Vitor, Eder T Rezende-Filho, Mariana F Crisóstomo, et al.
2021. “SRS: A Powerful r Package for Normalizing Whole-Genome Shotgun
Sequencing Data with Uneven Sequencing Depth.” *PeerJ* 9: e9593.
<https://doi.org/10.7717/peerj.9593>.

Legendre, Pierre, and Eugene D Gallagher. 2001. “Ecologically Meaningful
Transformations for Ordination of Species Data.” *Oecologia* 129 (2):
271–80. <https://doi.org/10.1007/s004420100716>.

Love, Michael I, Wolfgang Huber, and Simon Anders. 2014. “Moderated
Estimation of Fold Change and Dispersion for RNA-Seq Data with DESeq2.”
*Genome Biology* 15 (12): 550.
<https://doi.org/10.1186/s13059-014-0550-8>.

McKnight, Donald T, Roger Huerlimann, Deborah S Bower, Lin Schwarzkopf,
Ross A Alford, and Kyall R Zenger. 2019. “Methods for Normalizing
Microbiome Data: An Ecological Perspective.” *Methods in Ecology and
Evolution* 10 (3): 389–400. <https://doi.org/10.1111/2041-210X.13115>.

McMurdie, Paul J, and Susan Holmes. 2014. “Waste Not, Want Not: Why
Rarefying Microbiome Data Is Inadmissible.” *PLoS Computational Biology*
10 (4): e1003531. <https://doi.org/10.1371/journal.pcbi.1003531>.

Mikryukov, Vladimir, Olesya Dulya, Alexander Zizka, Mohammad Bahram,
Niloufar Hagh-Doust, Sten Anslan, Oleh Prylutskyi, et al. 2023.
“Connecting the Multiple Dimensions of Global Soil Fungal Diversity.”
*Science Advances* 9 (48): eadj8016.
<https://doi.org/10.1126/sciadv.adj8016>.

Paulson, Joseph N, O Colin Stine, Héctor Corrada Bravo, and Mihai Pop.
2013. “Differential Abundance Analysis for Microbial Marker-Gene
Surveys.” *Nature Methods* 10 (12): 1200–1202.
<https://doi.org/10.1038/nmeth.2658>.

Quinn, Thomas P, Ionas Erb, Mark F Richardson, and Tamsyn M Crowley.
2018. “Understanding Sequencing Data as Compositions: An Outlook and
Review.” *Bioinformatics* 34 (16): 2870–78.
<https://doi.org/10.1093/bioinformatics/bty175>.

Robinson, Mark D, and Alicia Oshlack. 2010. “A Scaling Normalization
Method for Differential Expression Analysis of RNA-Seq Data.” *Genome
Biology* 11 (3): R25. <https://doi.org/10.1186/gb-2010-11-3-r25>.

Weiss, Sophie, Zhenjiang Zech Xu, Shyamal Peddada, Amnon Amir, Kyle
Bittinger, Antonio Gonzalez, Catherine Lozupone, et al. 2017.
“Normalization and Microbial Differential Abundance Strategies Depend
Upon Data Characteristics.” *Microbiome* 5 (1): 27.
<https://doi.org/10.1186/s40168-017-0237-y>.
