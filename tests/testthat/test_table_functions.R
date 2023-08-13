data(data_fungi)

test_that("tax_datatable function works fine with data_fungi dataset", {
  expect_silent(tax_datatable(data_fungi))
  expect_silent(tax_datatable(data_fungi, taxonomic_level = c(1:2)))
  expect_s3_class(tax_datatable(data_fungi),"datatables")
  expect_s3_class(tax_datatable(data_fungi, taxonomic_level = c(1:2)),"datatables")
  expect_message(tax_datatable(data_fungi, modality = data_fungi@sam_data$Height))
  expect_s3_class(tax_datatable(data_fungi, modality = data_fungi@sam_data$Height),"datatables")
})

