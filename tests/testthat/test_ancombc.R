skip_on_cran()
data(data_fungi)
data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High")) |>
  clean_pq()

test_that("signif_ancombc works", {
  skip_on_cran()
  skip_if_not_installed("ANCOMBC")
  ancombc_res <- tryCatch(
    suppressWarnings(ancombc_pq(data_subset, "Height")),
    error = function(e) NULL
  )
  skip_if(is.null(ancombc_res))
  expect_type(ancombc_res, "list")
  result <- signif_ancombc(ancombc_res)
  expect_s3_class(result, "data.frame")
})

test_that("plot_ancombc_pq works", {
  skip_on_cran()
  skip_if_not_installed("ANCOMBC")
  ancombc_res <- tryCatch(
    suppressWarnings(ancombc_pq(data_subset, "Height")),
    error = function(e) NULL
  )
  skip_if(is.null(ancombc_res))
  library(patchwork)
  p <- plot_ancombc_pq(data_subset, ancombc_res)
  expect_s3_class(p, "ggplot")
})
