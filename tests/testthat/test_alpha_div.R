data(data_fungi)
data(enterotype)

test_that("hill_curves_pq works", {
  skip_on_cran()
  if (requireNamespace("iNEXT", quietly = TRUE)) {
    data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
    result <- hill_curves_pq(data_subset, "Height")
    expect_s3_class(result, "ggplot")
  }
})

test_that("hill_test_rarperm_pq works", {
  skip_on_cran()
  data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
  result <- hill_test_rarperm_pq(data_subset, "Height", nperm = 9)
  expect_type(result, "list")
})
