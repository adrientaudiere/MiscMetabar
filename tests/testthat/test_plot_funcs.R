data(data_fungi)
data(enterotype)

test_that("plot_complexity_pq works", {
  skip_on_cran()
  p <- plot_complexity_pq(data_fungi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_refseq_pq works", {
  skip_on_cran()
  p <- plot_refseq_pq(data_fungi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_refseq_extremity_pq works", {
  skip_on_cran()
  p <- plot_refseq_extremity_pq(data_fungi)
  expect_s3_class(p, "ggplot")
})

test_that("plot_seq_ratio_pq works", {
  skip_on_cran()
  p <- plot_seq_ratio_pq(data_fungi, "Height")
  expect_s3_class(p, "ggplot")
})

test_that("plot_var_part_pq works", {
  skip_on_cran()
  data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
  p <- plot_var_part_pq(data_subset, vec_variables = c("Height", "Time"))
  expect_s3_class(p, "ggplot")
})

test_that("accu_plot_balanced_modality works", {
  skip_on_cran()
  data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
  p <- accu_plot_balanced_modality(data_subset, "Height", nperm = 9, step = 1000)
  expect_s3_class(p, "ggplot")
})
