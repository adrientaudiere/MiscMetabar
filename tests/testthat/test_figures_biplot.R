data("data_fungi")
data("GlobalPatterns")
GP <- GlobalPatterns

data_fungi_2trees <- subset_samples(data_fungi, data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X"))
data_fungi_abun <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)

test_that("biplot_pq works", {
  expect_message(biplot_pq(data_fungi_2trees, merge_sample_by = "Tree_name"))
  expect_message(biplot_pq(data_fungi_2trees, fact = "Tree_name", merge_sample_by = "Tree_name"))
  expect_error(biplot_pq(data_fungi, merge_sample_by = "Tree_name"), "biplot_pq needs only two samples")
  expect_error(biplot_pq(data_fungi_2trees, fact = "Tree_name"), "biplot_pq needs only two samples")
  expect_error(biplot_pq(data_fungi_2trees, merge_sample_by = "tRREE_name"))
})


test_that("multi_biplot_pq works with data_fungi dataset", {
  p1 <- multi_biplot_pq(data_fungi_abun, split_by = "Time", na_remove = FALSE)
  p2 <- multi_biplot_pq(data_fungi_abun, "Height")
  expect_s3_class(p1[[1]], "ggplot")
  expect_type(p1, "list")
  expect_s3_class(p2[[1]], "ggplot")
  expect_type(p2, "list")
})