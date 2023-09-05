data("data_fungi")
data("data_fungi_sp_known")
data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns
data_fungi_2trees <- subset_samples(data_fungi, data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X"))
GP_archae <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
data_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")

test_that("hill_pq works with data_fungi dataset", {
  expect_message(expect_message(hill_pq(data_fungi, "Height")))
  expect_message(expect_message(hill_pq(data_fungi, "Height", add_points = TRUE)))
  expect_silent(suppressMessages(hill_pq(clean_pq(subset_samples_pq(data_fungi, !is.na(data_fungi@sam_data$Height))), "Height", letters = TRUE)))
  expect_silent(suppressMessages(hill_pq(data_fungi, "Height", add_points = TRUE, color_fac = "Time")))
  expect_equal(length(hill_pq(data_fungi, "Height", add_points = TRUE)), 4)
  expect_s3_class(hill_pq(data_fungi, "Height", add_points = TRUE)[[1]], "ggplot")
})

test_that("hill_pq works with GP dataset", {
  expect_message(hill_pq(GP, "SampleType"))
  expect_message(hill_pq(GP, "SampleType", add_points = TRUE))
  expect_silent(suppressMessages(hill_pq(GP, "SampleType", letters = TRUE)))
  expect_silent(suppressMessages(hill_pq(GP, "SampleType", add_points = TRUE)))
  expect_equal(length(hill_pq(GP, "SampleType", add_points = TRUE)), 4)
  expect_s3_class(hill_pq(GP, "SampleType", add_points = TRUE)[[1]], "ggplot")
})

test_that("accu_plot works with GlobalPatterns dataset", {
  expect_silent(accu_plot(GP_archae, fact = "X.SampleID", by.fact = TRUE))
  expect_silent(accu_plot(GP_archae, fact = "X.SampleID", by.fact = TRUE, print_sam_names = TRUE))
  expect_silent(accu_plot(GP_archae, "SampleType", add_nb_seq = TRUE, by.fact = TRUE))
  expect_silent(suppressWarnings(accu_plot(GP_archae, "SampleType", add_nb_seq = FALSE, by.fact = TRUE)))
  expect_error(accu_plot(GP_archae))
})

test_that("accu_plot works with data_fungi dataset", {
  expect_silent(accu_plot(data_basidio, fact = "Height", by.fact = TRUE))
  expect_error(accu_plot(data_basidio, fact = "Height", by.fact = FALSE))
  expect_silent(accu_plot(data_basidio, fact = "Height", by.fact = TRUE, print_sam_names = TRUE))
  expect_silent(accu_plot(data_basidio, "Height", add_nb_seq = TRUE, by.fact = TRUE))
  expect_silent(accu_plot(data_basidio, "Height", add_nb_seq = FALSE, by.fact = TRUE))
  expect_error(accu_plot(data_basidio))
})