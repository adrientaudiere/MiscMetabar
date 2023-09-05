
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
#> dplyr (1.1.2 -> 1.1.3) [CRAN]
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/tmp/RtmpluOtzG/remotes62e682f7ffc53/adrientaudiere-MiscMetabar-2cb7839/DESCRIPTION’ ... OK
#> * preparing ‘MiscMetabar’:
#> * checking DESCRIPTION meta-information ... OK
#> Avis : read_pq.Rd:26: unknown macro '\t'
#> Avis : write_pq.Rd:74: unknown macro '\t'
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘MiscMetabar_0.34.tar.gz’
```

You can install the developement version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("adrientaudiere/MiscMetabar", ref = "dev")
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/tmp/RtmpluOtzG/remotes62e682b7bdb84/adrientaudiere-MiscMetabar-1dc2e0c/DESCRIPTION’ ... OK
#> * preparing ‘MiscMetabar’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘MiscMetabar_0.40.tar.gz’
#> Warning in i.p(...): l'installation du package
#> '/tmp/RtmpluOtzG/file62e683fb2fa28/MiscMetabar_0.40.tar.gz' a eu un statut de
#> sortie non nul
```

## Some use of MiscMetabar

### Summarize a physeq object

``` r
library("MiscMetabar")
library("phyloseq")
library("magrittr")
data("data_fungi")
summary_plot_pq(data_fungi)
```

<img src="man/figures/README-example-1.png" width="100%" />

### Circle for visualize distribution of taxa in function of samples variables

``` r
circle_pq(data_fungi, "Height", taxa = "Class", add_nb_seq = F)
#> Only 10 taxa are plot (41.67%). Use 'min_prop_tax' to plot more taxa
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
