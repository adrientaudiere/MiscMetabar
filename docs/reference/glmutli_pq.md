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
  library("divent")
  res_glmulti <-
    glmutli_pq(data_fungi_mini,
      "Hill_0 ~ Hill_1 + Abundance + Time + Height",
      level = 1
    )
  res_glmulti
  res_glmulti_interaction <-
    glmutli_pq(data_fungi_mini,
      "Hill_0 ~ Abundance + Time + Height",
      level = 2
    )
  res_glmulti_interaction
}
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Joining with `by = join_by(Sample)`
#> Initialization...
#> TASK: Exhaustive screening of candidate set.
#> Fitting...
#> Completed.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Joining with `by = join_by(Sample)`
#> Initialization...
#> TASK: Exhaustive screening of candidate set.
#> Fitting...
#> 
#> After 50 models:
#> Best model: Hill_0~1+Abundance+Time+Time:Abundance+Height:Abundance
#> Crit= 380.266307886255
#> Mean crit= 439.149731147241
#> Completed.
#>                            estimates unconditional_interval nb_model
#> HeightHigh:Time         1.021021e-03           4.896053e-06        8
#> Abundance:HeightHigh    1.444574e-06           1.376218e-11        8
#> HeightLow              -2.818502e-01           2.975283e-01       32
#> HeightMiddle           -2.284000e-01           2.804621e-01       32
#> HeightLow:Time         -4.802362e-02           7.037018e-03       32
#> HeightMiddle:Time      -5.451488e-02           8.777216e-03       32
#> Abundance:HeightLow     9.952410e-05           2.196739e-08       32
#> Abundance:HeightMiddle  1.942462e-04           3.617651e-08       32
#> Time                    1.927954e-01           1.152254e-02       32
#> Abundance               4.367041e-04           3.508660e-08       32
#> Abundance:Time         -3.131657e-05           1.596068e-10       32
#>                         importance        alpha               variable
#> HeightHigh:Time        0.006596239 4.348738e-03        HeightHigh:Time
#> Abundance:HeightHigh   0.014891518 7.317402e-06   Abundance:HeightHigh
#> HeightLow              0.293275977 1.076681e+00              HeightLow
#> HeightMiddle           0.293275977 1.046605e+00           HeightMiddle
#> HeightLow:Time         0.409654935 1.655170e-01         HeightLow:Time
#> HeightMiddle:Time      0.409654935 1.848108e-01      HeightMiddle:Time
#> Abundance:HeightLow    0.671985280 2.935384e-04    Abundance:HeightLow
#> Abundance:HeightMiddle 0.671985280 3.755994e-04 Abundance:HeightMiddle
#> Time                   0.900282159 2.125915e-01                   Time
#> Abundance              0.933935731 3.708996e-04              Abundance
#> Abundance:Time         0.967098115 2.508155e-05         Abundance:Time
# }
```
