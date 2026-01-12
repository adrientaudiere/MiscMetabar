data(data_fungi)

data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High")) |>
  rarefy_even_depth() |> 
  clean_pq()

test_that("hill_curves_pq works", {
  if (requireNamespace("iNEXT", quietly = TRUE)) {
    result <- hill_curves_pq(data_subset, "Height")
    expect_s3_class(result, "ggplot")
  }
})

test_that("hill_test_rarperm_pq works", {
  result <- hill_test_rarperm_pq(data_subset, "Height", nperm = 9)
  expect_type(result, "list")
  expect_s3_class(result[[1]], "statsExpressions")
})
