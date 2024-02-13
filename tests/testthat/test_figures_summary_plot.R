data("data_fungi")
data("data_fungi_sp_known")
data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns

test_that("summary_plot_pq works with data_fungi dataset", {
  expect_message(summary_plot_pq(data_fungi))
  skip_on_cran()
  expect_s3_class(summary_plot_pq(data_fungi), "ggplot")
  expect_message(summary_plot_pq(data_fungi, add_info = FALSE))
  expect_message(summary_plot_pq(data_fungi, add_info = FALSE, min_seq_samples = 33))
  expect_message(summary_plot_pq(data_fungi) + scale_fill_viridis_d())
  expect_silent(summary_plot_pq(data_fungi, clean_pq = FALSE))
})

test_that("summary_plot_pq works with GP dataset", {
  expect_message(summary_plot_pq(GP))
  skip_on_cran()
  expect_s3_class(summary_plot_pq(GP), "ggplot")
  expect_message(summary_plot_pq(GP, add_info = FALSE))
  expect_message(summary_plot_pq(GP, add_info = FALSE, min_seq_samples = 33))
  expect_message(summary_plot_pq(GP) + scale_fill_viridis_d())
  expect_silent(summary_plot_pq(GP, clean_pq = FALSE))
})

test_that("summary_plot_pq works with enterotype dataset", {
  expect_message(summary_plot_pq(enterotype))
  skip_on_cran()
  expect_s3_class(summary_plot_pq(enterotype), "ggplot")
  expect_message(summary_plot_pq(enterotype, add_info = FALSE))
  expect_message(summary_plot_pq(enterotype, add_info = FALSE, min_seq_samples = 33))
  expect_message(summary_plot_pq(enterotype) + scale_fill_viridis_d())
  expect_silent(summary_plot_pq(enterotype, clean_pq = FALSE))
})
