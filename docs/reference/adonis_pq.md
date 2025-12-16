# Permanova on a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html)
function in the case of `physeq` object.

## Usage

``` r
adonis_pq(
  physeq,
  formula,
  dist_method = "bray",
  merge_sample_by = NULL,
  na_remove = FALSE,
  correction_for_sample_size = FALSE,
  rarefy_nb_seqs = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- formula:

  (required) the right part of a formula for
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html).
  Variables must be present in the `physeq@sam_data` slot.

- dist_method:

  (default "bray") the distance used. See
  [`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
  for all available distances or run
  [`phyloseq::distanceMethodList()`](https://rdrr.io/pkg/phyloseq/man/distanceMethodList.html).
  For aitchison and robust.aitchison distance,
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function is directly used.

- merge_sample_by:

  a vector to determine which samples to merge using the
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  function. Need to be in `physeq@sam_data`

- na_remove:

  (logical, default FALSE) If set to TRUE, remove samples with NA in the
  variables set in formula.

- correction_for_sample_size:

  (logical, default FALSE) If set to TRUE, the sample size (number of
  sequences by samples) is added to formula in the form
  `y~Library_Size + Biological_Effect` following recommendation of
  [Weiss et al.
  2017](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y).
  `correction_for_sample_size` overcome `rarefy_nb_seqs` if both are
  TRUE.

- rarefy_nb_seqs:

  (logical, default FALSE) Rarefy each sample (before merging if
  merge_sample_by is set) using
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html).
  if `correction_for_sample_size` is TRUE, rarefy_nb_seqs will have no
  effect.

- verbose:

  (logical, default TRUE) If TRUE, prompt some messages.

- ...:

  Additional arguments passed on to
  [`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html)
  function. Note that the parameter `by` is important. If by is set to
  NULL (default) the p-value is computed for the entire model. by = NULL
  will assess the overall significance of all terms together, by =
  "terms" will assess significance for each term (sequentially from
  first to last), setting by = "margin" will assess the marginal effects
  of the terms (each marginal term analyzed in a model with all other
  variables), by = "onedf" will analyze one-degree-of-freedom contrasts
  sequentially. The argument is passed on to anova.cca.

## Value

The function returns an anova.cca result object with a new column for
partial R^2. See help of
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/vegan-defunct.html)
for more information.

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
if you use this function.

## Author

Adrien Taudière

## Examples

``` r
data(enterotype)
# \donttest{
adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE)
#> Taxa are now in columns.
#> SeqTech
#> Enterotype
#> 9 were discarded due to NA in variables present in formula.
#> Permutation test for adonis under reduced model
#> Permutation: free
#> Number of permutations: 999
#> 
#> vegan::adonis2(formula = .formula, data = metadata)
#>           Df SumOfSqs      R2      F Pr(>F)    
#> Model      8   38.194 0.70766 79.278  0.001 ***
#> Residual 262   15.778 0.29234                  
#> Total    270   53.972 1.00000                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE, by = "terms")
#> Taxa are now in columns.
#> SeqTech
#> Enterotype
#> 9 were discarded due to NA in variables present in formula.
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> vegan::adonis2(formula = .formula, data = metadata, by = "terms")
#>                     Df SumOfSqs      R2        F Pr(>F)    
#> SeqTech              2   29.175 0.54055 242.2289  0.001 ***
#> Enterotype           2    8.651 0.16028  71.8250  0.001 ***
#> SeqTech:Enterotype   4    0.368 0.00683   1.5293  0.096 .  
#> Residual           262   15.778 0.29234                    
#> Total              270   53.972 1.00000                    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE, by = "onedf")
#> Taxa are now in columns.
#> SeqTech
#> Enterotype
#> 9 were discarded due to NA in variables present in formula.
#> Permutation test for adonis under reduced model
#> Sequential test for contrasts
#> Permutation: free
#> Number of permutations: 999
#> 
#> vegan::adonis2(formula = .formula, data = metadata, by = "onedf")
#>                             Df SumOfSqs       R2        F Pr(>F)    
#> SeqTechPyro454               1   28.796  0.53353 478.1659  0.001 ***
#> SeqTechSanger                1    0.379  0.00702   6.2918  0.002 ** 
#> Enterotype2                  1    4.282  0.07933  71.1013  0.001 ***
#> Enterotype3                  1    4.369  0.08095  72.5487  0.001 ***
#> SeqTechPyro454:Enterotype2   1   -0.556 -0.01030  -9.2315  1.000    
#> SeqTechSanger:Enterotype2    1    0.232  0.00430   3.8560  0.011 *  
#> SeqTechPyro454:Enterotype3   1    0.443  0.00820   7.3535  0.001 ***
#> SeqTechSanger:Enterotype3    1    0.249  0.00462   4.1394  0.009 ** 
#> Residual                   262   15.778  0.29234                    
#> Total                      270   53.972  1.00000                    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE, by = "margin")
#> Taxa are now in columns.
#> SeqTech
#> Enterotype
#> 9 were discarded due to NA in variables present in formula.
#> Permutation test for adonis under reduced model
#> Marginal effects of terms
#> Permutation: free
#> Number of permutations: 999
#> 
#> vegan::adonis2(formula = .formula, data = metadata, by = "margin")
#>                     Df SumOfSqs      R2      F Pr(>F)
#> SeqTech:Enterotype   4    0.368 0.00683 1.5293  0.123
#> Residual           262   15.778 0.29234              
#> Total              270   53.972 1.00000              

adonis_pq(enterotype, "SeqTech", dist_method = "jaccard", by = "terms")
#> Taxa are now in columns.
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> vegan::adonis2(formula = .formula, data = metadata, by = "terms")
#>           Df SumOfSqs      R2      F Pr(>F)    
#> SeqTech    2   31.330 0.40211 93.147  0.001 ***
#> Residual 277   46.585 0.59789                  
#> Total    279   77.915 1.00000                  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
adonis_pq(enterotype, "SeqTech", dist_method = "robust.aitchison", by = "terms")
#> Taxa are now in columns.
#> Permutation test for adonis under reduced model
#> Terms added sequentially (first to last)
#> Permutation: free
#> Number of permutations: 999
#> 
#> vegan::adonis2(formula = .formula, data = metadata, by = "terms")
#>           Df SumOfSqs      R2    F Pr(>F)    
#> SeqTech    2   121403 0.95598 3008  0.001 ***
#> Residual 277     5590 0.04402                
#> Total    279   126992 1.00000                
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# }
```
