# Partition the Variation of a phyloseq object by 2, 3, or 4 Explanatory Matrices

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

The function partitions the variation in otu_table using distance (Bray
per default) with respect to two, three, or four explanatory tables,
using adjusted R² in redundancy analysis ordination (RDA) or
distance-based redundancy analysis. If response is a single vector,
partitioning is by partial regression. Collinear variables in the
explanatory tables do NOT have to be removed prior to partitioning. See
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html)
for more information.

## Usage

``` r
var_par_pq(
  physeq,
  list_component,
  dist_method = "bray",
  dbrda_computation = TRUE
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- list_component:

  (required) A named list of 2, 3 or four vectors with names from the
  `@sam_data` slot.

- dist_method:

  (default "bray") the distance used. See
  [`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
  for all available distances or run
  [`phyloseq::distanceMethodList()`](https://rdrr.io/pkg/phyloseq/man/distanceMethodList.html).
  For "aitchison" and "robust.aitchison" distance,
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function is directly used.

- dbrda_computation:

  (logical) Do dbrda computations are runned for each individual
  component (each name of the list component) ?

## Value

an object of class "varpart", see
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html)

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html)
if you use this function.

## See also

[`var_par_rarperm_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_rarperm_pq.md),
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html),
[`plot_var_part_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_var_part_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("vegan")) {
  data_fungi_woNA <-
    subset_samples(data_fungi, !is.na(Time) & !is.na(Height))
  res_var <- var_par_pq(data_fungi_woNA,
    list_component = list(
      "Time" = c("Time"),
      "Size" = c("Height", "Diameter")
    ),
    dbrda_computation = TRUE
  )
}
# }
```
