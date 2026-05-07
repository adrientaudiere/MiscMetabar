# Rules

## Rules for the package development

### Documentation

- Always indicate required params `(required)`
- Indicate default values only when (i) a set of values is possible
  \[e.g `"bray"` for the parameter `method` in
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)\]
  or (ii) when this value is well thought out and a good default value
  for most users \[e.g. the number of permutation or processors, or the
  level of identity to cluster at 97%\].
- Indicate the type for *logical* and *integer* params.
- Homogenize the params names across function.
- Prefer ASVs to OTUs denomination, even thought both are almost
  interchangeable in the MiscMetabar package.

### Lifecycle

- **Experimental**: First status for a function
- **Maturing**: Some tests and/or analyses make stronger these functions
- **Stable**: Good level of confidence

### Tests and examples

#### Special cases of external softwares

Some examples and test required the installation of external softwares.
I used MiscMetabar internal functions (such as
\[is_vsearch_installed()\]) to conditionnaly test
(`if(is_vsearch_installed()){...test...}`) and run examples
`#' @examplesIf is_vsearch_installed()`.

#### Special case of CRAN

I use `\donttest{}` for some long examples (in roxygen documentation)
and
[`testthat::skip_on_cran()`](https://testthat.r-lib.org/reference/skip.html)
in long test.

#### Special case of windows

Some tests and examples are not tested on windows. I used
`testthat::skip_on_os("windows")` inside test and
`#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"` in
examples (roxygen documentation). Here is a list of functions with some
limitations or not working at all on windows OS:

- [`build_phytree_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/build_phytree_pq.md)
- [`count_seq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/count_seq.md)
- [`cutadapt_remove_primers()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/cutadapt_remove_primers.md)
- [`krona()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/krona.md)
- [`merge_krona()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/merge_krona.md)
- [`multipatt_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/multipatt_pq.md)
- [`plot_tsne_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/plot_tsne_pq.md)
- [`rotl_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/rotl_pq.md)
- [`save_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/save_pq.md)
- [`tax_datatable()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/tax_datatable.md)
- [`track_wkflow()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/track_wkflow.md)
- [`track_wkflow_samples()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/track_wkflow_samples.md)
- [`tsne_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/tsne_pq.md)
- [`venn_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/venn_pq.md)

## Session information

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Pop!_OS 24.04 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Paris
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
    ##  [5] xfun_0.57         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
    ##  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
    ## [13] pkgdown_2.2.0     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
    ## [17] compiler_4.5.2    tools_4.5.2       ragg_1.5.2        bslib_0.10.0     
    ## [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
    ## [25] rlang_1.2.0       fs_2.0.1          htmlwidgets_1.6.4
