# Automated model selection and multimodel inference with (G)LMs for phyloseq

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

See [`glmulti::glmulti()`](https://rdrr.io/pkg/glmulti/man/glmulti.html)
for more information.

## Usage

``` r
glmutli_pq(
  physeq,
  formula,
  fitfunction = "lm",
  hill_scales = c(0, 1, 2),
  aic_step = 2,
  confsetsize = 100,
  plotty = FALSE,
  level = 1,
  method = "h",
  crit = "aicc",
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- formula:

  (required) a formula for
  [`glmulti::glmulti()`](https://rdrr.io/pkg/glmulti/man/glmulti.html)
  Variables must be present in the `physeq@sam_data` slot or be one of
  hill number defined in hill_scales or the variable Abundance which
  refer to the number of sequences per sample.

- fitfunction:

  (default "lm")

- hill_scales:

  (a vector of integer) The list of q values to compute the hill number
  H^q. If Null, no hill number are computed. Default value compute the
  Hill number 0 (Species richness), the Hill number 1 (exponential of
  Shannon Index) and the Hill number 2 (inverse of Simpson Index).

- aic_step:

  The value between AIC scores to cut for.

- confsetsize:

  The number of models to be looked for, i.e. the size of the returned
  confidence set.

- plotty:

  (logical) Whether to plot the progress of the IC profile when running.

- level:

  If 1, only main effects (terms of order 1) are used to build the
  candidate set. If 2, pairwise interactions are also used (higher order
  interactions are currently ignored)

- method:

  The method to be used to explore the candidate set of models. If "h"
  (default) an exhaustive screening is undertaken. If "g" the genetic
  algorithm is employed (recommended for large candidate sets). If "l",
  a very fast exhaustive branch-and-bound algorithm is used. Package
  leaps must then be loaded, and this can only be applied to linear
  models with covariates and no interactions. If "d", a simple summary
  of the candidate set is printed, including the number of candidate
  models.

- crit:

  The Information Criterion to be used. Default is the small-sample
  corrected AIC (aicc). This should be a function that accepts a fitted
  model as first argument. Other provided functions are the classic AIC,
  the Bayes IC (bic), and QAIC/QAICc (qaic and qaicc).

- ...:

  Additional arguments passed on to
  [`glmulti::glmulti()`](https://rdrr.io/pkg/glmulti/man/glmulti.html)
  function

## Value

A data.frame summarizing the glmulti results with columns

-estimates -unconditional_interval -nb_model" -importance -alpha

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`glmulti::glmulti()`](https://rdrr.io/pkg/glmulti/man/glmulti.html) if
you use this function.

## See also

[`glmulti::glmulti()`](https://rdrr.io/pkg/glmulti/man/glmulti.html)

## Examples

``` r
# \donttest{
if (requireNamespace("glmulti")) {
  res_glmulti <-
    glmutli_pq(data_fungi, "Hill_0 ~ Hill_1 + Abundance + Time + Height", level = 1)
  res_glmulti
  res_glmulti_interaction <-
    glmutli_pq(data_fungi, "Hill_0 ~ Abundance + Time + Height", level = 2)
  res_glmulti
}
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Initialization...
#> TASK: Exhaustive screening of candidate set.
#> Fitting...
#> Completed.
#> Taxa are now in rows.
#> Joining with `by = join_by(Sample)`
#> Initialization...
#> TASK: Exhaustive screening of candidate set.
#> Fitting...
#> 
#> After 50 models:
#> Best model: Hill_0~1+Abundance+Time+Time:Abundance+Height:Abundance+Height:Time
#> Crit= 1069.11608982306
#> Mean crit= 1218.19009955263
#> Completed.
#>                estimates unconditional_interval nb_model importance
#> Hill_1       3.062117997           1.868174e-01        8          1
#> Abundance    0.002959644           8.478374e-08        8          1
#> Time         0.789091999           2.443263e-01        8          1
#> HeightLow    6.884340946           3.444196e+01        8          1
#> HeightMiddle 0.339123798           3.727962e+01        8          1
#>                     alpha     variable
#> Hill_1       8.570200e-01       Hill_1
#> Abundance    5.773492e-04    Abundance
#> Time         9.800932e-01         Time
#> HeightLow    1.163660e+01    HeightLow
#> HeightMiddle 1.210648e+01 HeightMiddle
# }
```
