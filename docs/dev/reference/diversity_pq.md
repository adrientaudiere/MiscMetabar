# Compute alpha- or beta-diversity metrics on a phyloseq object via ecodive

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wrapper around
[`ecodive::beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.html)
and
[`ecodive::alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.html)
providing fast, C-accelerated alpha- and beta-diversity calculations
directly on a `phyloseq` object, plus support for `"robust_aitchison"`
(alias `"robust.aitchison"`).

## Usage

``` r
diversity_pq(
  physeq,
  metric = "bray",
  norm = "none",
  pseudocount = NULL,
  impute = FALSE,
  add_to_phyloseq = FALSE,
  nproc = 1L,
  ...
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- metric:

  (character, default `"bray"`) Name(s) of the diversity metric(s) to
  compute. See
  [`ecodive::list_metrics()`](https://cmmr.github.io/ecodive/reference/metrics.html)
  for all supported alpha and beta metric IDs, plus `"robust_aitchison"`
  (or `"robust.aitchison"`). All elements of `metric` must belong to the
  same family (either all `"beta"` or all `"alpha"`).

- norm:

  (character, default `"none"`) Normalization passed to
  [`ecodive::beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.html)
  (`"none"`, `"percent"`, `"binary"`, `"clr"`). For alpha-diversity
  metrics that require proportions, `"none"` is mapped to
  [`ecodive::alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.html)'s
  default `"percent"`.

- pseudocount:

  (numeric or `NULL`, default `NULL`) Pseudocount added to counts to
  handle zeros when `metric = "aitchison"` or `norm = "clr"`.

- impute:

  (logical, default `FALSE`) Only used when `metric` is
  `"robust_aitchison"` or `"robust.aitchison"`. If `FALSE` (default),
  uses `ecodive`'s zero-filled rCLR Euclidean distance. If `TRUE`, uses
  `vegan::vegdist(..., method = "robust.aitchison")` with OptSpace
  low-rank matrix completion.

- add_to_phyloseq:

  (logical, default `FALSE`) Only used for alpha-diversity metrics. If
  `TRUE`, appends the computed alpha-diversity columns to
  `physeq@sam_data` and returns the updated `phyloseq` object.

- nproc:

  (integer, default `1L`) Number of CPU threads passed to `ecodive`
  (`cpus`). Defaults to `1L` to comply with CRAN check limits.

- ...:

  Additional arguments forwarded to
  [`ecodive::beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.html),
  [`ecodive::alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.html),
  or
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  (when `impute = TRUE`).

## Value

- For **beta-diversity**: a
  [stats::dist](https://rdrr.io/r/stats/dist.html) object when
  `length(metric) == 1`, or a named list of
  [stats::dist](https://rdrr.io/r/stats/dist.html) objects when
  `length(metric) > 1`.

- For **alpha-diversity**: a
  [tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
  with column `Sample` and one column per metric (when
  `add_to_phyloseq = FALSE`), or a `phyloseq` object with those columns
  added to `sample_data(physeq)` (when `add_to_phyloseq = TRUE`).

## Details

**Robust Aitchison distance (`metric = "robust_aitchison"`):** Two
definitions exist in the microbiome literature (both citing Martino et
al. 2019):

- **`impute = FALSE` (default, `ecodive` definition):** applies the
  robust centered log-ratio (rCLR) transformation to strictly positive
  counts in each sample (`log(x) - mean(log(x[x > 0]))` for `x > 0`, and
  `0` otherwise) and computes the Euclidean distance via `ecodive`. This
  avoids pseudocounts and runs in milliseconds in C, without iterative
  matrix completion.

- **`impute = TRUE` (DEICODE / `vegan` OptSpace definition):** applies
  rCLR leaving zeros as `NA`, completes the low-rank matrix with
  OptSpace
  ([`vegan::decostand()`](https://vegandevs.github.io/vegan/reference/decostand.html)
  with `method = "rclr", impute = TRUE`, called by
  `vegan::vegdist(..., method = "robust.aitchison")`), and computes the
  Euclidean distance on the completed matrix. Because `vegan` already
  implements OptSpace matrix completion natively, no external `rpca`
  dependency is needed.

## References

Martino C, Morton JT, Marotz CA, Thompson LR, Tripathi A, Knight R,
Zengler K (2019). A novel sparse compositional technique reveals
microbial perturbations. *mSystems* 4(1): e00016-19.
[doi:10.1128/mSystems.00016-19](https://doi.org/10.1128/mSystems.00016-19)

## See also

[`ecodive::beta_div()`](https://cmmr.github.io/ecodive/reference/beta_div.html),
[`ecodive::alpha_div()`](https://cmmr.github.io/ecodive/reference/alpha_div.html),
[`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html),
[`hill_pq()`](https://adrientaudiere.github.io/MiscMetabar/dev/reference/hill_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("ecodive", quietly = TRUE)) {
  pq_5 <- prune_samples(sample_names(data_fungi_mini)[1:5], data_fungi_mini)
  # Beta-diversity (Bray-Curtis and Robust Aitchison)
  diversity_pq(pq_5, metric = "bray")
  diversity_pq(pq_5, metric = "robust_aitchison")

  # Alpha-diversity
  diversity_pq(pq_5, metric = c("shannon", "simpson", "observed"))
}
#> Taxa are now in columns.
#> # A tibble: 5 × 4
#>   Sample                         shannon simpson observed
#>   <chr>                            <dbl>   <dbl>    <dbl>
#> 1 A10-005-B_S188_MERGED.fastq.gz  1.33   0.699          6
#> 2 A10-005-H_S189_MERGED.fastq.gz  0.0118 0.00274        5
#> 3 A10-005-M_S190_MERGED.fastq.gz  1.42   0.727          7
#> 4 A12-007_S191_MERGED.fastq.gz    0.506  0.185         10
#> 5 A12-007-B_S2_MERGED.fastq.gz    0.687  0.494          2
# }
```
