data(data_fungi)

test_that("tax_datatable function works fine with data_fungi dataset", {
  expect_silent(tax_datatable(data_fungi))
  expect_silent(tax_datatable(data_fungi, taxonomic_level = c(1:2)))
  expect_s3_class(tax_datatable(data_fungi),"datatables")
  expect_s3_class(tax_datatable(data_fungi, taxonomic_level = c(1:2)),"datatables")
  expect_warning(tax_datatable(data_fungi, modality = data_fungi@sam_data$Height))
  expect_s3_class(tax_datatable(data_fungi, modality = data_fungi@sam_data$Height),"datatables")
})

data(enterotype)

test_that("tax_datatable function works fine with enterotype dataset", {
  expect_silent(tax_datatable(enterotype))
  expect_s3_class(tax_datatable(enterotype),"datatables")
  expect_silent(tax_datatable(enterotype, modality = enterotype@sam_data$SeqTech))
  expect_s3_class(tax_datatable(enterotype, modality = enterotype@sam_data$SeqTech),"datatables")
})

