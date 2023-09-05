data(data_fungi)
data("enterotype")

test_that("dist_bycol works fine", {
  expect_equal(length(dist_bycol(data_fungi@otu_table, as_binary_otu_table(data_fungi)@otu_table)), 2)
  expect_error(length(dist_bycol(data_fungi@otu_table, enterotype@otu_table)))
})

test_that("all_object_size works fine", {
  expect_type(all_object_size(), "double")
})
