# Test multiple times effect of factor on Hill diversity with different rarefaction even depth

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

This reduce the risk of a random drawing of a exceptional situation of
an unique rarefaction.

## Usage

``` r
hill_test_rarperm_pq(
  physeq,
  fact,
  hill_scales = c(0, 1, 2),
  nperm = 99,
  sample.size = min(sample_sums(physeq)),
  verbose = FALSE,
  progress_bar = TRUE,
  p_val_signif = 0.05,
  type = "non-parametrique",
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor in `physeq@sam_data` used to plot
  different lines

- hill_scales:

  (a vector of integer) The list of q values to compute the hill number
  H^q. If Null, no hill number are computed. Default value compute the
  Hill number 0 (Species richness), the Hill number 1 (exponential of
  Shannon Index) and the Hill number 2 (inverse of Simpson Index).

- nperm:

  (int) The number of permutations to perform.

- sample.size:

  (int) A single integer value equal to the number of reads being
  simulated, also known as the depth. See
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html).

- verbose:

  (logical). If TRUE, print additional information.

- progress_bar:

  (logical, default TRUE) Do we print progress during the calculation?

- p_val_signif:

  (float, `[0:1]`) The mimimum value of p-value to count a test as
  significant int the `prop_signif` result.

- type:

  A character specifying the type of statistical approach (See
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html)
  for more details):

  - "parametric"

  - "nonparametric"

  - "robust"

  - "bayes"

- ...:

  Additional arguments passed on to
  [`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html)
  function

## Value

A list of 6 components :

- method

- expressions

- plots

- pvals

- prop_signif

- statistics

## See also

[`ggstatsplot::ggbetweenstats()`](https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html),
[`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("ggstatsplot")) {
  hill_test_rarperm_pq(data_fungi, "Time", nperm = 2)
  res <- hill_test_rarperm_pq(data_fungi, "Height", nperm = 9, p.val = 0.9)
  patchwork::wrap_plots(res$plots[[1]])
  res$plots[[1]][[1]] + res$plots[[2]][[1]] + res$plots[[3]][[1]]
  res$prop_signif
  res_para <- hill_test_rarperm_pq(data_fungi, "Height", nperm = 9, type = "parametrique")
  res_para$plots[[1]][[1]] + res_para$plots[[2]][[1]] + res_para$plots[[3]][[1]]
  res_para$pvals
  res_para$method
  res_para$expressions[[1]]
}
#> 
  |                                                        
  |                                                  |   0%
  |                                                        
  |=========================                         |  50%
  |                                                        
  |==================================================| 100%
  |                                                        
  |                                                  |   0%
  |                                                        
  |======                                            |  11%
  |                                                        
  |===========                                       |  22%
  |                                                        
  |=================                                 |  33%
  |                                                        
  |======================                            |  44%
  |                                                        
  |============================                      |  56%
  |                                                        
  |=================================                 |  67%
  |                                                        
  |=======================================           |  78%
  |                                                        
  |============================================      |  89%
  |                                                        
  |==================================================| 100%
  |                                                        
  |                                                  |   0%
  |                                                        
  |======                                            |  11%
  |                                                        
  |===========                                       |  22%
  |                                                        
  |=================                                 |  33%
  |                                                        
  |======================                            |  44%
  |                                                        
  |============================                      |  56%
  |                                                        
  |=================================                 |  67%
  |                                                        
  |=======================================           |  78%
  |                                                        
  |============================================      |  89%
  |                                                        
  |==================================================| 100%
#> list(italic("F")["Welch"](2, 84.92) == "0.08", italic(p) == "0.92", 
#>     widehat(omega["p"]^2) == "0.00", CI["95%"] ~ "[" * "0.00", 
#>     "1.00" * "]", italic("n")["obs"] == "131")
# }
```
