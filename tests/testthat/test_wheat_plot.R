test_that("wheat_plot returns a ggplot", {
  set.seed(42)
  df <- data.frame(value = rnorm(200, mean = 50, sd = 10))
  p <- suppressWarnings(wheat_plot(df, value, binwidth = 2))
  expect_s3_class(p, "ggplot")
})

test_that("wheat_plot chooses a binwidth automatically", {
  set.seed(1)
  df <- data.frame(value = rnorm(100))
  p <- suppressWarnings(wheat_plot(df, value))
  expect_s3_class(p, "ggplot")
  # one point per observation
  expect_equal(nrow(p$data), nrow(df))
})

test_that("wheat_plot works on taxa_sums of a phyloseq object", {
  df <- data.frame(value = taxa_sums(data_fungi_mini))
  p <- suppressWarnings(wheat_plot(df, value, binwidth = 2000))
  expect_s3_class(p, "ggplot")
})

test_that("wheat_plot errors on a non-numeric column", {
  df <- data.frame(value = letters[1:10])
  expect_error(suppressWarnings(wheat_plot(df, value)), "numeric")
})
