# Resolve conflict in a vector of taxonomy values

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

Internally used in the function
[`assign_blastn()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_blastn.md)
with method="vote" and
[`assign_vsearch_lca()`](https://adrientaudiere.github.io/MiscMetabar/reference/assign_vsearch_lca.md)
if `top_hits_only` is FALSE and `vote_algorithm` is not NULL.

## Usage

``` r
resolve_vector_ranks(
  vec,
  method = c("consensus", "rel_majority", "abs_majority", "preference", "unanimity"),
  strict = FALSE,
  second_method = c("consensus", "rel_majority", "abs_majority", "unanimity"),
  nb_agree_threshold = 1,
  preference_index = NULL,
  collapse_string = "/",
  replace_collapsed_rank_by_NA = FALSE
)
```

## Arguments

- vec:

  (required) A vector of (taxonomic) values

- method:

  One of "consensus", "rel_majority", "abs_majority", "preference" or
  "unanimity". See details.

- strict:

  (logical, default FALSE). If TRUE, NA are considered as informative in
  resolving conflict (i.e. NA are taking into account in vote). See
  details for more informations.

- second_method:

  One of "consensus", "rel_majority", "abs_majority", or "unanimity".
  Only used if method = "preference". See details.

- nb_agree_threshold:

  (Int, default 1) The minimum number of times a value must arise to be
  selected using vote. If 2, we only kept taxonomic value present at
  least 2 times in the vector.

- preference_index:

  (Int. default NULL). Required if method="preference". Useless for
  other method. The preference index is the index of the value in vec
  for which we have a preference.

- collapse_string:

  (default '/'). The character to collapse taxonomic names when multiple
  assignment is done.

- replace_collapsed_rank_by_NA:

  (logical, default FALSE). If set to TRUE, all multiple assignments
  (all taxonomic rank including the 'collapse_string' parameter) are
  replaced by NA.

## Value

a vector of length 1 (one character value)

## Details

- `unanimity`: Only keep taxonomic value when all methods are agree

  - By default, the value with NA are not taking into account
    (strict=FALSE)

  - If `strict` , one NA in the row is sufficient to return a NA

- `consensus`: Keep all taxonomic values separated by a '/' (separation
  can be modify using param `collapse_string`)

  - If `strict` is TRUE, NA are also added to taxonomic vector such as
    'Tiger/Cat/NA' instead of 'Tiger/Cat'

- `abs_majority`: Keep the most found taxonomic value if it represent at
  least half of all taxonomic values.

  - If `strict` is TRUE, NA values are also count to determine the
    majority. For example, a vector of taxonomic rank c("A", "A", "A",
    "B", NA, NA) will give a value of 'A' if `strict` is FALSE (default)
    but a value of NA if `strict` is TRUE.

- `rel_majority`: Keep the most found taxonomic value. In case of
  equality, apply a consensus strategy (collapse values separated by a
  '/') across the most found taxonomic values.

  - If `strict` is TRUE, NA are considered as a rank and can win the
    relative majority vote. If `strict` is FALSE (default), NA are
    removed before ranking the taxonomic values.

  - `nb_agree_threshold`: Only keep return value when at least x methods
    agreed with x is set by parameter `nb_agree_threshold`. By default,
    (`nb_agree_threshold` = 1): a majority of one is enough.

- `preference`: Keep the value from a preferred column.

  - when the value is NA in the preferred column, apply a second
    strategy (by default `consensus`) to the other column (parameter
    `second_method`). Note that the parameters `strict` and
    `nb_agree_threshold` are used for the second_method consensus.

## Author

Adrien Taudière

## Examples

``` r
resolve_vector_ranks(c("A"))
#> [1] "A"
resolve_vector_ranks(c("A"),
  method = "preference",
  preference_index = 1
)
#> [1] "A"
resolve_vector_ranks(c("A"), method = "abs_majority")
#> [1] "A"
resolve_vector_ranks(c("A"), method = "rel_majority")
#> [1] "A"
resolve_vector_ranks(c("A"),
  method = "rel_majority",
  nb_agree_threshold = 2
)
#> [1] NA
resolve_vector_ranks(c("A"), method = "unanimity")
#> [1] "A"

resolve_vector_ranks(c("A", "A", "A"))
#> [1] "A"
resolve_vector_ranks(c("A", "A", "A"),
  method = "preference",
  preference_index = 1
)
#> [1] "A"
resolve_vector_ranks(c("A", "A", "A"), method = "abs_majority")
#> [1] "A"
resolve_vector_ranks(c("A", "A", "A"), method = "rel_majority")
#> [1] "A"
resolve_vector_ranks(c("A", "A", "A"), method = "unanimity")
#> [1] "A"

resolve_vector_ranks(c(NA, NA, NA))
#> [1] NA
resolve_vector_ranks(c(NA, NA, NA),
  method = "preference",
  preference_index = 1
)
#> [1] NA
resolve_vector_ranks(c(NA, NA, NA), method = "abs_majority")
#> [1] NA
resolve_vector_ranks(c(NA, NA, NA), method = "rel_majority")
#> [1] NA
resolve_vector_ranks(c(NA, NA, NA), method = "unanimity")
#> [1] NA

resolve_vector_ranks(c("A", "A", NA))
#> [1] "A"
resolve_vector_ranks(c("A", "A", NA),
  method = "preference",
  preference_index = 1
)
#> [1] "A"
resolve_vector_ranks(c("A", "A", NA), method = "rel_majority")
#> [1] "A"
resolve_vector_ranks(c("A", "A", NA), method = "abs_majority")
#> [1] "A"
resolve_vector_ranks(c("A", "A", NA, NA),
  method = "abs_majority",
  strict = FALSE
)
#> [1] "A"
resolve_vector_ranks(c("A", "A", NA, NA),
  method = "abs_majority",
  strict = TRUE
)
#> [1] NA
resolve_vector_ranks(c("A", "A", NA), method = "unanimity")
#> [1] "A"
resolve_vector_ranks(c("A", "A", NA),
  method = "unanimity",
  strict = TRUE
)
#> [1] NA

resolve_vector_ranks(c("A", "B", NA))
#> [1] "A/B"
resolve_vector_ranks(c("A", "B", NA), strict = TRUE)
#> [1] "A/B/NA"
resolve_vector_ranks(c("A", "B", NA),
  method = "preference",
  preference_index = 1
)
#> [1] "A"
resolve_vector_ranks(c("A", "B", NA), method = "abs_majority")
#> [1] NA
resolve_vector_ranks(c("A", "B", NA), method = "rel_majority")
#> [1] "A/B"
resolve_vector_ranks(c("A", "B", NA),
  method = "rel_majority",
  strict = TRUE
)
#> [1] "A/B/NA"
resolve_vector_ranks(c("A", "B", NA),
  method = "rel_majority",
  nb_agree_threshold = 2
)
#> [1] "A/B"
resolve_vector_ranks(c("A", "B", NA), method = "unanimity")
#> [1] NA

resolve_vector_ranks(c("A", NA, NA))
#> [1] "A"
resolve_vector_ranks(c("A", NA, NA), method = "rel_majority")
#> [1] "A"
resolve_vector_ranks(c("A", NA, NA), method = "unanimity")
#> [1] "A"
resolve_vector_ranks(c("A", NA, NA),
  method = "preference",
  preference_index = 1
)
#> [1] "A"
resolve_vector_ranks(c("A", NA, NA),
  method = "preference",
  preference_index = 2
)
#> [1] "A"
resolve_vector_ranks(c("A", NA, "B"),
  method = "preference",
  preference_index = 2
)
#> [1] "A/B"
resolve_vector_ranks(c("A", NA, "B"),
  method = "preference",
  preference_index = 2, second_method = "abs_majority"
)
#> [1] NA

resolve_vector_ranks(c("A", "B", "B"))
#> [1] "A/B"
resolve_vector_ranks(c("A", "B", "B"),
  method = "preference",
  preference_index = 1
)
#> [1] "A"
resolve_vector_ranks(c("A", "B", "B"), method = "abs_majority")
#> [1] "B"
resolve_vector_ranks(c("A", "B", "B"), method = "rel_majority")
#> [1] "B"
resolve_vector_ranks(c("A", "B", "B"), method = "unanimity")
#> [1] NA


resolve_vector_ranks(c("A", "A", "A", "B", NA, NA))
#> [1] "A/B"
resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
  strict = TRUE
)
#> [1] "A/B/NA"
resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
  method = "abs_majority"
)
#> [1] "A"
resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
  method = "abs_majority",
  strict = TRUE
)
#> [1] NA
resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
  method = "preference", preference_index = 6, second_method = "abs_majority"
)
#> [1] "A"
resolve_vector_ranks(c("A", "A", "A", "B", NA, NA, NA),
  method = "preference", preference_index = 6, second_method = "abs_majority"
)
#> [1] "A"
resolve_vector_ranks(c("A", "A", "A", "B", NA, NA, NA),
  method = "preference", preference_index = 6, second_method = "abs_majority",
  strict = TRUE
)
#> [1] NA
```
