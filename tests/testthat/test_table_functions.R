data(data_fungi)

test_that("tax_datatable function works fine with data_fungi dataset", {
  expect_silent(taxdt <- tax_datatable(data_fungi))
  expect_s3_class(taxdt, "datatables")
  expect_silent(taxdt <- tax_datatable(data_fungi, taxonomic_level = c(1:2)))
  expect_s3_class(taxdt, "datatables")
  expect_no_warning(taxdt <- suppressWarnings(tax_datatable(data_fungi, modality = data_fungi@sam_data$Height)))
  expect_s3_class(taxdt, "datatables")
})

data(enterotype)

test_that("tax_datatable function works fine with enterotype dataset", {
  expect_silent(tax_datatable(enterotype))
  expect_s3_class(tax_datatable(enterotype), "datatables")
  expect_silent(tax_datatable(enterotype, modality = enterotype@sam_data$SeqTech))
  expect_s3_class(tax_datatable(enterotype, modality = enterotype@sam_data$SeqTech), "datatables")
})


data_fungi_low_high <- subset_samples(data_fungi, Height %in% c("Low", "High"))
data_fungi_low_high_withNA <- data_fungi_low_high
data_fungi_low_high_withNA@sam_data[["Height"]][1] <- NA

test_that(" compare_pairs_pq function works fine with data_fungi dataset", {
  expect_s3_class(compare_pairs_pq(data_fungi_low_high, bifactor = "Height", merge_sample_by = "Height"), "tbl_df")
  expect_message(expect_message(compare_pairs_pq(data_fungi_low_high_withNA, bifactor = "Height", merge_sample_by = "Height")))
  expect_equal(dim(compare_pairs_pq(data_fungi_low_high, bifactor = "Height", merge_sample_by = "Height")), c(1, 13))
  expect_s3_class(compare_pairs_pq(data_fungi_low_high, bifactor = "Height", merge_sample_by = "Height", nb_min_seq = 2), "tbl_df")
  expect_s3_class(compare_pairs_pq(data_fungi_low_high, bifactor = "Height", merge_sample_by = "Height", veg_index = "simpson"), "tbl_df")
  expect_s3_class(compare_pairs_pq(data_fungi_low_high_withNA, bifactor = "Height", merge_sample_by = "Height", modality = "Time"), "tbl_df")
  expect_equal(dim(compare_pairs_pq(data_fungi_low_high_withNA, bifactor = "Height", merge_sample_by = "Height", modality = "Time")), c(4, 13))
})
