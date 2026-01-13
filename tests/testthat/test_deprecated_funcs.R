data(data_fungi)

test_that("clean_physeq deprecated function works", {
  expect_warning(result <- clean_physeq(data_fungi), "deprecated")
  expect_s4_class(result, "phyloseq")
})
