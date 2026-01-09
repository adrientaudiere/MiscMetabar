data(data_fungi)
data(enterotype)

test_that("signif_ancombc works", {
  skip_on_cran()
  if (requireNamespace("ANCOMBC", quietly = TRUE)) {
    data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
    # Run ancombc first
    ancombc_res <- ancombc_pq(data_subset, "Height")
    result <- signif_ancombc(ancombc_res)
    expect_s3_class(result, "data.frame")
  }
})

test_that("plot_ancombc_pq works", {
  skip_on_cran()
  if (requireNamespace("ANCOMBC", quietly = TRUE)) {
    data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
    ancombc_res <- ancombc_pq(data_subset, "Height")
    p <- plot_ancombc_pq(ancombc_res, data_subset)
    expect_s3_class(p, "ggplot")
  }
})
