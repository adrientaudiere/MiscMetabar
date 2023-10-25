
![R](https://img.shields.io/badge/r-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
<a href="https://zenodo.org/badge/latestdoi/268765075"><img src="https://zenodo.org/badge/268765075.svg" alt="DOI"></a>
[![codecov](https://codecov.io/gh/adrientaudiere/MiscMetabar/graph/badge.svg?token=NXFRSIKYC0)](https://codecov.io/gh/adrientaudiere/MiscMetabar)
[![Contributor
Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](code_of_conduct.md)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- devtools::build_readme() -->

# MiscMetabar

See the pkdown site
[here](https://adrientaudiere.github.io/MiscMetabar/).

Biological studies, especially in ecology, health sciences and taxonomy,
need to describe the biological composition of samples. During the last
twenty years, (i) the development of DNA sequencing, (ii) reference
databases, (iii) high-throughput sequencing (HTS), and (iv)
bioinformatics resources have allowed the description of biological
communities through metabarcoding. Metabarcoding involves the sequencing
of millions (*meta*-) of short regions of specific DNA (*-barcoding*,
Valentini, Pompanon, and Taberlet (2009)) often from environmental
samples (eDNA, Taberlet et al. (2012)) such as human stomach contents,
lake water, soil and air.

`MiscMetabar` aims to facilitate the **description**,
**transformation**, **exploration** and **reproducibility** of
metabarcoding analysis using R. The development of `MiscMetabar` relies
heavily on the R packages
[`dada2`](https://benjjneb.github.io/dada2/index.html),
[`phyloseq`](https://joey711.github.io/phyloseq/) and
[`targets`](https://books.ropensci.org/targets/).

## Installation

There is no CRAN or bioconductor version of MiscMetabar for now (work in
progress).

You can install the stable version from [GitHub](https://github.com/)
with:

``` r
install.packages("devtools")
devtools::install_github("adrientaudiere/MiscMetabar")
```

You can install the developement version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("adrientaudiere/MiscMetabar", ref = "dev")
```

## Some use of MiscMetabar

See vignettes in the
[MiscMetabar](https://adrientaudiere.github.io/MiscMetabar/) website for
more examples.

### Summarize a physeq object

``` r
library("MiscMetabar")
library("phyloseq")
library("magrittr")
data("data_fungi")
summary_plot_pq(data_fungi)
```

<img src="man/figures/README-example-1.png" width="100%" />

### Alpha-diversity analysis

``` r
p <- MiscMetabar::hill_pq(data_fungi, variable = "Height")
#> Taxa are now in rows.
#> Cleaning suppress 0 taxa and 0 samples.
p$plot_Hill_0
```

<div class="figure">

<img src="man/figures/README-unnamed-chunk-4-1.png" alt="Hill number 1" width="100%" />
<p class="caption">
Hill number 1
</p>

</div>

``` r
p$plot_tuckey
```

<div class="figure">

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="Result of the Tuckey post-hoc test" width="100%" />
<p class="caption">
Result of the Tuckey post-hoc test
</p>

</div>

### Beta-diversity analysis

``` r
ggvenn_pq(data_fungi, fact = "Height") +
  ggplot2::scale_fill_distiller(palette = "BuPu", direction = 1) +
  labs(title = "Share number of ASV among Height in tree")
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

### Installation of other softwares for debian Linux distributions

#### blastn

``` sh
sudo apt-get install ncbi-blast+
```

#### vsearch

``` sh
sudo apt-get install vsearch
```

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-taberlet2012" class="csl-entry">

Taberlet, Pierre, Eric Coissac, Mehrdad Hajibabaei, and Loren H
Rieseberg. 2012. “Environmental Dna.” *Molecular Ecology*. Wiley Online
Library. <https://doi.org/10.1002/(issn)2637-4943>.

</div>

<div id="ref-valentini2009" class="csl-entry">

Valentini, Alice, François Pompanon, and Pierre Taberlet. 2009. “DNA
Barcoding for Ecologists.” *Trends in Ecology & Evolution* 24 (2):
110–17. <https://doi.org/10.1016/j.tree.2008.09.011>.

</div>

</div>
