test_that("resolve_vector_ranks works with consensus method", {
  expect_equal(resolve_vector_ranks(c("A")), "A")
  expect_equal(resolve_vector_ranks(c("A", "A", "A")), "A")
  expect_equal(resolve_vector_ranks(c("A", "B")), "A/B")
  expect_equal(resolve_vector_ranks(c("A", "B", NA)), "A/B")
  expect_equal(resolve_vector_ranks(c("A", "B", NA), strict = TRUE), "A/B/NA")
  expect_true(is.na(resolve_vector_ranks(c(NA, NA, NA))))
})

test_that("resolve_vector_ranks works with rel_majority method", {
  expect_equal(resolve_vector_ranks(c("A"), method = "rel_majority"), "A")
  expect_equal(resolve_vector_ranks(c("A", "A", "A"), method = "rel_majority"), "A")
  expect_equal(resolve_vector_ranks(c("A", "B", "B"), method = "rel_majority"), "B")
  expect_equal(resolve_vector_ranks(c("A", "B", NA), method = "rel_majority"), "A/B")
  expect_equal(resolve_vector_ranks(c("A", "B", NA), method = "rel_majority", strict = TRUE), "NA")
  skip_on_cran()
  expect_equal(resolve_vector_ranks(c("A", "B", NA), method = "rel_majority", nb_agree_threshold = 2), NA)
  expect_equal(resolve_vector_ranks(c("A", NA, NA), method = "rel_majority"), "A")
})

test_that("resolve_vector_ranks works with abs_majority method", {
  expect_equal(resolve_vector_ranks(c("A"), method = "abs_majority"), "A")
  expect_equal(resolve_vector_ranks(c("A", "A", "A"), method = "abs_majority"), "A")
  expect_equal(resolve_vector_ranks(c("A", "B", "B"), method = "abs_majority"), "B")
  expect_true(is.na(resolve_vector_ranks(c("A", "B", NA), method = "abs_majority")))
  skip_on_cran()
  expect_equal(resolve_vector_ranks(c("A", "A", NA), method = "abs_majority"), "A")
  expect_true(is.na(resolve_vector_ranks(c("A", "A", NA, NA), method = "abs_majority", strict = TRUE)))
  expect_equal(resolve_vector_ranks(c("A", "A", NA, NA), method = "abs_majority", strict = FALSE), "A")
})

test_that("resolve_vector_ranks works with unanimity method", {
  expect_equal(resolve_vector_ranks(c("A"), method = "unanimity"), "A")
  expect_equal(resolve_vector_ranks(c("A", "A", "A"), method = "unanimity"), "A")
  expect_true(is.na(resolve_vector_ranks(c("A", "B", "B"), method = "unanimity")))
  expect_true(is.na(resolve_vector_ranks(c("A", "B", NA), method = "unanimity")))
  skip_on_cran()
  expect_equal(resolve_vector_ranks(c("A", "A", NA), method = "unanimity"), "A")
  expect_true(is.na(resolve_vector_ranks(c("A", "A", NA), method = "unanimity", strict = TRUE)))
})

test_that("resolve_vector_ranks works with preference method", {
  expect_equal(resolve_vector_ranks(c("A"), method = "preference", preference_index = 1), "A")
  expect_equal(
    resolve_vector_ranks(c("A", "B", "B"), method = "preference", preference_index = 1),
    "A"
  )
  skip_on_cran()
  expect_equal(
    resolve_vector_ranks(c("A", NA, NA), method = "preference", preference_index = 2),
    "A"
  )
  expect_equal(
    resolve_vector_ranks(c("A", NA, "B"), method = "preference", preference_index = 2),
    "A/B"
  )
  expect_equal(
    resolve_vector_ranks(c("A", NA, "B"),
      method = "preference",
      preference_index = 2, second_method = "abs_majority"
    ),
    NA
  )
  expect_error(resolve_vector_ranks(c("A"), method = "preference"))
})

test_that("resolve_vector_ranks works with collapse_string parameter", {
  expect_equal(
    resolve_vector_ranks(c("A", "B"), collapse_string = "-"),
    "A-B"
  )
  skip_on_cran()
  expect_equal(
    resolve_vector_ranks(c("A", "B"), collapse_string = "; "),
    "A; B"
  )
})

test_that("resolve_vector_ranks works with replace_collapsed_rank_by_NA", {
  expect_true(is.na(
    resolve_vector_ranks(c("A", "B"), replace_collapsed_rank_by_NA = TRUE)
  ))
  expect_equal(
    resolve_vector_ranks(c("A", "A"), replace_collapsed_rank_by_NA = TRUE),
    "A"
  )
})
