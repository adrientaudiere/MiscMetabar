data(data_fungi)
data(enterotype)

test_that("adonis_rarperm_pq works", {
  skip_on_cran()
  data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
  result <- adonis_rarperm_pq(data_subset, "Height", nperm = 9)
  expect_type(result, "list")
})

test_that("var_par_pq works", {
  skip_on_cran()
  data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
  result <- var_par_pq(data_subset, vec_variables = c("Height", "Time"))
  expect_type(result, "list")
})

test_that("var_par_rarperm_pq works", {
  skip_on_cran()
  data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
  result <- var_par_rarperm_pq(data_subset, vec_variables = c("Height", "Time"), nperm = 9)
  expect_type(result, "list")
})
