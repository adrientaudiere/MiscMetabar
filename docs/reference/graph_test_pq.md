# Performs graph-based permutation tests on phyloseq object

[![lifecycle-maturing](https://img.shields.io/badge/lifecycle-maturing-blue)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper of
[`phyloseqGraphTest::graph_perm_test()`](https://rdrr.io/pkg/phyloseqGraphTest/man/graph_perm_test.html)
for quick plot with important statistics

## Usage

``` r
graph_test_pq(
  physeq,
  fact,
  merge_sample_by = NULL,
  nperm = 999,
  return_plot = TRUE,
  title = "Graph Test",
  na_remove = FALSE,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (required) Name of the factor to cluster samples by modalities. Need
  to be in `physeq@sam_data`. This should be a factor with two or more
  levels.

- merge_sample_by:

  a vector to determine which samples to merge using
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  function. Need to be in `physeq@sam_data`

- nperm:

  (int) The number of permutations to perform.

- return_plot:

  (logical) Do we return only the result of the test, or do we plot the
  result?

- title:

  The title of the Graph.

- na_remove:

  (logical, default FALSE) If set to TRUE, remove samples with NA in the
  variables set in formula.

- ...:

  Other params for be passed on to
  [`phyloseqGraphTest::graph_perm_test()`](https://rdrr.io/pkg/phyloseqGraphTest/man/graph_perm_test.html)
  function

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)2 plot
with a subtitle indicating the pvalue and the number of permutations

## Details

This function is mainly a wrapper of the work of others. Please cite
`phyloseqGraphTest` package.

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("phyloseqGraphTest")) {
  data(enterotype)
  graph_test_pq(enterotype, fact = "SeqTech")
  graph_test_pq(enterotype, fact = "Enterotype", na_remove = TRUE)
}
#> 9 were discarded due to NA in variables present in formula.

# }
```
