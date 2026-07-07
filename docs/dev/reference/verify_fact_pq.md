# Verify that grouping columns exist in the `sam_data` slot of a phyloseq object

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Check that the column name(s) supplied as a grouping factor (`fact` or
its synonym `modality`) or as a two-level grouping factor (`bifactor`)
are actually present in the
[`sample_data()`](https://rdrr.io/pkg/phyloseq/man/sample_data-methods.html)
(`sam_data`) slot of a `phyloseq` object. For `bifactor`, the column is
additionally required to have **exactly two** levels.

The goal is to raise a clear, early error listing the available columns
instead of letting a missing column propagate into a cryptic downstream
failure (e.g. `NULL` being coerced to a zero-level factor). Most
functions of the pqverse that take a `physeq` together with a `fact`,
`modality` or `bifactor` argument call this helper internally.

Each argument accepts a character vector, so several columns can be
checked in a single call. `NULL` arguments are skipped.

## Usage

``` r
verify_fact_pq(
  physeq,
  fact = NULL,
  bifactor = NULL,
  modality = NULL,
  call = rlang::caller_env()
)
```

## Arguments

- physeq:

  (required) a
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- fact:

  (default NULL) Name(s) of the `sam_data` column(s) used as a primary
  grouping factor. Checked for presence only.

- bifactor:

  (default NULL) Name(s) of the `sam_data` column(s) used as a two-level
  grouping factor. Checked for presence **and** for having exactly two
  levels.

- modality:

  (default NULL) Synonym of `fact` kept for the functions of the pqverse
  that name their grouping argument `modality`. Checked for presence
  only.

- call:

  (default
  [`rlang::caller_env()`](https://rlang.r-lib.org/reference/stack.html))
  The calling environment, used to point the error message at the
  user-facing function rather than at `verify_fact_pq()` itself.

## Value

Invisibly returns `physeq`. The function is called for its side effect:
it throws an informative error when a check fails.

## Author

Adrien Taudière

## Examples

``` r
# Presence check (passes silently)
verify_fact_pq(data_fungi_mini, fact = "Height")

if (FALSE) { # \dontrun{
# Missing column: error lists the available sam_data columns
verify_fact_pq(data_fungi_mini, fact = "Heigth")

# bifactor with more than two levels: error (Height has 3 levels)
verify_fact_pq(data_fungi_mini, bifactor = "Height")
} # }

# A genuine two-level column passes the bifactor check
data_2h <- subset_samples_pq(
  data_fungi_mini,
  data_fungi_mini@sam_data$Height %in% c("Low", "High")
)
verify_fact_pq(data_2h, bifactor = "Height")
```
