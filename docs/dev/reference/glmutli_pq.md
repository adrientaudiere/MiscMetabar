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
  q = c(0, 1, 2),
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
  hill number defined in q or the variable Abundance which refer to the
  number of sequences per sample.

- fitfunction:

  (default "lm")

- q:

  (a vector of integer) The list of q values to compute the hill number
  H^q. If Null, no hill number are computed. Default value compute the
  Hill number 0 (Species richness), the Hill number 1 (exponential of
  Shannon Index) and the Hill number 2 (inverse of Simpson Index). Hill
  numbers are more appropriate in DNA metabarcoding studies when `q > 0`
  (Alberdi & Gilbert, 2019; Calderón-Sanou et al., 2019).

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

  (character, default aicc) The Information Criterion to be used.
  Default is the small-sample corrected AIC (aicc). This should be a
  function that accepts a fitted model as first argument. Other provided
  functions are the classic AIC, the Bayes IC (bic), and QAIC/QAICc
  (qaic and qaicc).

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

## References

Alberdi, A., & Gilbert, M. T. P. (2019). A guide to the application of
Hill numbers to DNA-based diversity analyses. *Molecular Ecology
Resources*.
[doi:10.1111/1755-0998.13014](https://doi.org/10.1111/1755-0998.13014)

Calderón-Sanou, I., Münkemüller, T., Boyer, F., Zinger, L., & Thuiller,
W. (2019). From environmental DNA sequences to ecological conclusions:
How strong is the influence of methodological choices? *Journal of
Biogeography*, 47.
[doi:10.1111/jbi.13681](https://doi.org/10.1111/jbi.13681)

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
  res_glmulti_interaction
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
#> Best model: Hill_0~1+Abundance+Time+Time:Abundance+Height:Abundance
#> Crit= 1162.46935121017
#> Mean crit= 1326.57756179615
#> Completed.
#>                            estimates unconditional_interval nb_model importance
#> HeightHigh:Time         0.0238428433           5.369345e-03        8 0.02339703
#> Abundance:HeightHigh    0.0001567326           9.040694e-08        8 0.05936981
#> HeightLow               1.3672307044           2.837356e+01       32 0.28630897
#> HeightMiddle           -2.8823946528           4.305119e+01       32 0.28630897
#> HeightLow:Time          0.6976197936           1.991032e+00       32 0.48361405
#> HeightMiddle:Time      -0.3490561261           1.442196e+00       32 0.48361405
#> Abundance:HeightLow     0.0009609861           1.351266e-06       32 0.52347724
#> Abundance:HeightMiddle  0.0008434925           1.380631e-06       32 0.52347724
#> Time                    1.9148085653           3.424966e+00       32 0.72978508
#> Abundance:Time         -0.0001497892           9.036339e-09       32 0.86061128
#> Abundance               0.0040684019           1.456572e-06       32 0.94026666
#>                               alpha               variable
#> HeightHigh:Time        1.442951e-01        HeightHigh:Time
#> Abundance:HeightHigh   5.895794e-04   Abundance:HeightHigh
#> HeightLow              1.050820e+01              HeightLow
#> HeightMiddle           1.292959e+01           HeightMiddle
#> HeightLow:Time         2.778215e+00         HeightLow:Time
#> HeightMiddle:Time      2.367383e+00      HeightMiddle:Time
#> Abundance:HeightLow    2.286482e-03    Abundance:HeightLow
#> Abundance:HeightMiddle 2.312397e-03 Abundance:HeightMiddle
#> Time                   3.645902e+00                   Time
#> Abundance:Time         1.874648e-04         Abundance:Time
#> Abundance              2.378267e-03              Abundance
# }
```
