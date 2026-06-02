# Test and plot multipatt result

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

A wrapper for the
[`indicspecies::multipatt()`](https://emf-creaf.github.io/indicspecies/reference/multipatt.html)
function in the case of `physeq` object.

## Usage

``` r
multipatt_pq(
  physeq,
  fact,
  p_adjust_method = "BH",
  pval = 0.05,
  control = permute::how(nperm = 999),
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

- p_adjust_method:

  (chr, default "BH"): the method used to adjust p-value

- pval:

  (int, default 0.05): the value to determine the significance of LCBD

- control:

  see `?indicspecies::multipatt()`

- ...:

  Additional arguments passed on to
  [`indicspecies::multipatt()`](https://emf-creaf.github.io/indicspecies/reference/multipatt.html)
  function

## Value

A ggplot2 object

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`indicspecies::multipatt()`](https://emf-creaf.github.io/indicspecies/reference/multipatt.html)
if you use this function.

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("indicspecies")) {
  multipatt_pq(subset_samples(data_fungi_mini, !is.na(Time)),
    fact = "Time", control = permute::how(nperm = 99))
  multipatt_pq(subset_samples(data_fungi_mini, !is.na(Time)),
    fact = "Time",
    max.order = 1, control = permute::how(nperm = 99)
  )
}
#> Taxa are now in columns.
#> Taxa are now in columns.

# }
```
