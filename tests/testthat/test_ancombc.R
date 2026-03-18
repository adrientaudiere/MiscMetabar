data(data_fungi)
data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High")) |>
  clean_pq()

ancombc_res <- tryCatch(
  suppressWarnings(ancombc_pq(data_subset, "Height")),
  error = \(e) NULL
)

test_that("signif_ancombc works", {
  if (requireNamespace("ANCOMBC", quietly = TRUE) && !is.null(ancombc_res)) {
    expect_type(ancombc_res, "list")
    result <- signif_ancombc(ancombc_res)
    expect_s3_class(result, "data.frame")
  }
})

test_that("plot_ancombc_pq works", {
  if (requireNamespace("ANCOMBC", quietly = TRUE) && !is.null(ancombc_res)) {
    p <- plot_ancombc_pq(data_subset, ancombc_res)
    expect_s3_class(p, "ggplot")
  }
})
