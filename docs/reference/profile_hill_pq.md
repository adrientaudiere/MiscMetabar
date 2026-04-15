# Hill diversity profile for a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Wraps
[`divent::profile_hill()`](https://ericmarcon.github.io/divent/reference/profile_hill.html)
to compute a Hill diversity profile across diversity orders for each
sample in a `phyloseq` object, and returns a ggplot2 object via
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html).

## Usage

``` r
profile_hill_pq(physeq, orders = seq(0, 2, 0.1), merge_sample_by = NULL, ...)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- orders:

  (numeric vector) Hill diversity orders to compute. Default
  `seq(0, 2, 0.1)`.

- merge_sample_by:

  (character or NULL) If not NULL, merge samples using
  [`merge_samples2()`](https://adrientaudiere.github.io/MiscMetabar/reference/merge_samples2.md)
  before computing profiles.

- ...:

  Additional arguments passed to
  [`divent::profile_hill()`](https://ericmarcon.github.io/divent/reference/profile_hill.html).

## Value

A ggplot2 object.

## See also

[`divent::profile_hill()`](https://ericmarcon.github.io/divent/reference/profile_hill.html),
[`hill_curves_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/hill_curves_pq.md)

## Examples

``` r
# \donttest{
profile_hill_pq(data_fungi_mini)
#> Taxa are now in columns.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> ! Sample coverage is 0, most estimators will return `NaN`.
#> Computing Hill profiles ■■■■■■■■■■■■■■■                   45% | ETA:  5s
#> Computing Hill profiles ■■■■■■■■■■■■■■■■■■■■■■■           73% | ETA:  4s

# }
```
