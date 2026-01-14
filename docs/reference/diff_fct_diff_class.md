# Compute different functions for different class of vector.

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Mainly an internal function useful in "sapply(..., tapply)" methods

## Usage

``` r
diff_fct_diff_class(
  x,
  numeric_fonction = mean,
  logical_method = "TRUE_if_one",
  character_method = "unique_or_na",
  ...
)
```

## Arguments

- x:

  : a vector

- numeric_fonction:

  : a function for numeric vector. For ex. `sum` or `mean`

- logical_method:

  : A method for logical vector. One of :

  - TRUE_if_one (default)

  - NA_if_not_all_TRUE

  - FALSE_if_not_all_TRUE

- character_method:

  : A method for character vector (and factor). One of :

  - unique_or_na (default)

  - more_frequent

  - more_frequent_without_equality

- ...:

  Additional arguments passed on to the numeric function (ex.
  na.rm=TRUE)

## Value

a single value

## Author

Adrien Taudière

## Examples

``` r
diff_fct_diff_class(
  data_fungi@sam_data$Sample_id,
  numeric_fonction = sum,
  na.rm = TRUE
)
#> [1] 17852
diff_fct_diff_class(
  data_fungi@sam_data$Time,
  numeric_fonction = mean,
  na.rm = TRUE
)
#> [1] 5.802469
diff_fct_diff_class(
  data_fungi@sam_data$Height == "Low",
  logical_method = "TRUE_if_one"
)
#> [1] TRUE
diff_fct_diff_class(
  data_fungi@sam_data$Height == "Low",
  logical_method = "NA_if_not_all_TRUE"
)
#> [1] NA
diff_fct_diff_class(
  data_fungi@sam_data$Height == "Low",
  logical_method = "FALSE_if_not_all_TRUE"
)
#> [1] FALSE
diff_fct_diff_class(
  data_fungi@sam_data$Height,
  character_method = "unique_or_na"
)
#> [1] NA
diff_fct_diff_class(
  c("IE", "IE"),
  character_method = "unique_or_na"
)
#> [1] "IE"
diff_fct_diff_class(
  c("IE", "IE", "TE", "TE"),
  character_method = "more_frequent"
)
#> [1] "IE"
diff_fct_diff_class(
  c("IE", "IE", "TE", "TE"),
  character_method = "more_frequent_without_equality"
)
#> [1] NA
```
