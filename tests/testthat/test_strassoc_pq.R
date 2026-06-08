data("data_fungi_mini", package = "MiscMetabar")
df_h <- subset_samples(data_fungi_mini, !is.na(Height))

test_that("strassoc_pq returns a tidy tibble without bootstrap", {
  skip_if_not_installed("indicspecies")
  res <- strassoc_pq(df_h, fact = "Height", func = "IndVal.g")
  expect_s3_class(res, "tbl_df")
  expect_true("taxon" %in% colnames(res))
  expect_identical(nrow(res), ntaxa(df_h))
})

test_that("strassoc_pq accepts other association indices", {
  skip_if_not_installed("indicspecies")
  expect_s3_class(strassoc_pq(df_h, fact = "Height", func = "A"), "tbl_df")
  expect_s3_class(strassoc_pq(df_h, fact = "Height", func = "B"), "tbl_df")
})

test_that("strassoc_pq returns a list of CI matrices with bootstrap", {
  skip_if_not_installed("indicspecies")
  res <- strassoc_pq(df_h, fact = "Height", func = "IndVal.g", nboot_ci = 99)
  expect_type(res, "list")
  expect_true(all(c("lowerCI", "stat", "upperCI") %in% names(res)))
})

test_that("strassoc_pq errors when the factor has fewer than two levels", {
  skip_if_not_installed("indicspecies")
  df_one <- subset_samples(df_h, Height == "Low")
  expect_error(strassoc_pq(df_one, fact = "Height"), "two levels")
})
