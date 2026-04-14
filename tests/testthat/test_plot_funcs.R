data(data_fungi)

test_that("plot_complexity_pq works", {
  p <- plot_complexity_pq(data_fungi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_refseq_pq works", {
  p <- plot_refseq_pq(data_fungi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_refseq_extremity_pq works", {
  p <- plot_refseq_extremity_pq(data_fungi)
  expect_s3_class(p[[1]], "ggplot")
  expect_s3_class(p[[2]], "ggplot")
  expect_length(p, 4)
})

test_that("plot_seq_ratio_pq works", {
  p <- plot_seq_ratio_pq(data_fungi)
  expect_s3_class(p, "ggplot")
})

test_that("accu_plot_balanced_modality works", {
  data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
  suppressWarnings(
    p <- accu_plot_balanced_modality(
      data_subset,
      "Height",
      nperm = 9,
      step = 10000
    )
  )
  expect_s3_class(p, "ggplot")
})
