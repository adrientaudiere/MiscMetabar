# Build a sample information tibble from physeq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Hill numbers are the number of equiprobable species giving the same
diversity value as the observed distribution.

Note that contrary to
[`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_pq.md),
this function does not take into account for difference in the number of
sequences per samples/modalities. You may use rarefy_by_sample = TRUE if
the mean number of sequences per samples differs among modalities.

## Usage

``` r
psmelt_samples_pq(
  physeq,
  q = c(0, 1, 2),
  hill_scales = lifecycle::deprecated(),
  filter_zero = TRUE,
  rarefy_by_sample = FALSE,
  rngseed = FALSE,
  verbose = TRUE,
  taxa_ranks = NULL,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- q:

  (numeric vector) Hill diversity orders to compute. If NULL, no Hill
  numbers are computed. Default computes Hill number 0 (species
  richness), 1 (exponential of Shannon index) and 2 (inverse of Simpson
  index). Formerly `q`. Hill numbers are more appropriate in DNA
  metabarcoding studies when `q > 0` (Alberdi & Gilbert, 2019;
  Calderón-Sanou et al., 2019).

- hill_scales:

  **\[deprecated\]** Use `q` instead.

- filter_zero:

  (logical, default TRUE) Do we filter non present OTU from samples ?
  For the moment, this has no effect on the result because the dataframe
  is grouped by samples with abundance summed across OTU.

- rarefy_by_sample:

  (logical, default FALSE) If TRUE, rarefy samples using
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  function.

- rngseed:

  (Optional). A single integer value passed to
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html),
  which is used to fix a seed for reproducibly random number generation
  (in this case, reproducibly random subsampling). If set to FALSE, then
  no fiddling with the RNG seed is performed, and it is up to the user
  to appropriately call set.seed beforehand to achieve reproducible
  results. Default is FALSE.

- verbose:

  (logical). If TRUE, print additional information.

- taxa_ranks:

  A vector of taxonomic ranks. For examples c("Family","Genus"). If taxa
  ranks is not set (default value = NULL), taxonomic information are not
  present in the resulting tibble.

- ...:

  Additional arguments passed to
  [`divent_hill_matrix_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/divent_hill_matrix_pq.md)
  and hence to
  [`divent::div_hill()`](https://ericmarcon.github.io/divent/reference/div_hill.html)
  (e.g. `estimator = "naive"`). Only used when `q` is not NULL.

## Value

A tibble with a row for each sample. Columns provide information from
`sam_data` slot as well as hill numbers, Abundance (nb of sequences),
and Abundance_log10 (*log10(1+Abundance)*).

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

## Author

Adrien Taudière

## Examples

``` r
if (requireNamespace("ggstatsplot")) {
  psm_tib <- psmelt_samples_pq(data_fungi_mini, hill_scales = c(0, 2, 7))
  ggstatsplot::ggbetweenstats(psm_tib, Height, Hill_0)
  ggstatsplot::ggbetweenstats(psm_tib, Height, Hill_7)

  psm_tib_tax <- psmelt_samples_pq(data_fungi_mini, taxa_ranks = c("Class", "Family"))
  ggplot(filter(psm_tib_tax, Abundance > 2000), aes(y = Family, x = Abundance, fill = Time)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Height)
}
#> Warning: The `hill_scales` argument of `psmelt_samples_pq()` is deprecated as of
#> MiscMetabar 0.15.1.
#> ℹ Please use the `q` argument instead.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Joining with `by = join_by(Sample)`
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Joining with `by = join_by(Sample)`
```
