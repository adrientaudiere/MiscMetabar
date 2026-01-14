data(data_fungi)
data(enterotype)

test_that("taxa_as_columns works", {
  result <- taxa_as_columns(data_fungi)
  expect_s4_class(result, "phyloseq")
  expect_false(taxa_are_rows(result))
})

test_that("taxa_as_rows works", {
  result <- taxa_as_rows(data_fungi)
  expect_s4_class(result, "phyloseq")
  expect_true(taxa_are_rows(result))
})

test_that("taxa_only_in_one_level works", {
  suppressWarnings(result <- taxa_only_in_one_level(data_fungi, "Height", "Low"))
  expect_equal(length(result), 124)
})

test_that("tbl_sum_taxtable works", {
  if (requireNamespace("gtsummary", quietly = TRUE)) {
    result <- tbl_sum_taxtable(data_fungi)
    expect_s3_class(result, "tbl_summary")
  }
})
