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
