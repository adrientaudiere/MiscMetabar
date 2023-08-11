data("data_fungi")
data("data_fungi_sp_known")
data("GlobalPatterns")

data_fungi_2trees <- subset_samples(data_fungi, data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X"))
GP_archae <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")

test_that("biplot_pq works", {
  expect_message(biplot_pq(data_fungi_2trees, merge_sample_by = "Tree_name"))
  expect_message(biplot_pq(data_fungi_2trees, fact = "Tree_name", merge_sample_by = "Tree_name"))
  expect_error(biplot_pq(data_fungi, merge_sample_by = "Tree_name"), "biplot_pq needs only two samples")
  expect_error(biplot_pq(data_fungi_2trees, fact = "Tree_name"), "biplot_pq needs only two samples")
  expect_error(biplot_pq(data_fungi_2trees, merge_sample_by = "tRREE_name"))
})

test_that("circle_pq works", {
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", nproc = 4))
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", min_prop_tax = 0.0001, min_prop_mod = 0.0001, nproc = 4))
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", nproc = 4, add_nb_seq = FALSE))
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", nproc = 4, rarefy = TRUE))
  expect_error(circle_pq(data_fungi_2trees, fact = "tRREE_name"))
})

test_that("graph_test_pq works", {
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name"))
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name", return_plot = FALSE))
  expect_message(graph_test_pq(subset_samples(data_fungi, !is.na(data_fungi@sam_data$Time)), fact = "Time", merge_sample_by = "Tree_name"))
  expect_error(graph_test_pq(data_fungi, fact = "Height"))
})

test_that("accu_plot works", {
  expect_silent(accu_plot(GP_archae, fact = "X.SampleID", by.fact = T))
  expect_silent(accu_plot(GP_archae, "SampleType", add_nb_seq = TRUE, by.fact = TRUE))
  expect_warning(accu_plot(GP_archae, "SampleType", add_nb_seq = FALSE, by.fact = TRUE))
  expect_error(accu_plot(GP_archae))
})
