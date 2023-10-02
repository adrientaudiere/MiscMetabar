
<a href="https://zenodo.org/badge/latestdoi/268765075"><img src="https://zenodo.org/badge/268765075.svg" alt="DOI"></a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- devtools::build_readme() -->

# MiscMetabar

The goal of MiscMetabar is to complete the great packages
[dada2](https://benjjneb.github.io/dada2/index.html),
[phyloseq](https://joey711.github.io/phyloseq/) and
[targets](https://books.ropensci.org/targets/). See the pkdown site
[here](https://adrientaudiere.github.io/MiscMetabar/).

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
