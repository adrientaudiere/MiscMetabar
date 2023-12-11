data("GlobalPatterns", package = "phyloseq")
GlobalPatterns@sam_data[, "Soil_logical"] <-
  ifelse(GlobalPatterns@sam_data[, "SampleType"] == "Soil", "Soil", "Not Soil")
test_that("hill_tuckey_pq function works fine with GlobalPatterns dataset", {
  expect_silent(suppressMessages(hill_tuckey_pq(GlobalPatterns, "Soil_logical")))
  expect_silent(suppressMessages(hill_tuckey_pq(GlobalPatterns, "SampleType")))
  expect_message(hill_tuckey_pq(GlobalPatterns, "SampleType", silent = FALSE))
  expect_s3_class(hill_tuckey_pq(GlobalPatterns, "SampleType"), "ggplot")
  expect_error(hill_tuckey_pq(GlobalPatterns, "SampleTYPE"))
})


data("enterotype")
test_that("hill_tuckey_pq function works fine with enterotype dataset", {
  expect_silent(hill_tuckey_pq(enterotype, "Nationality"))
  expect_message(hill_tuckey_pq(enterotype, "Nationality", silent = FALSE))
  expect_s3_class(hill_tuckey_pq(enterotype, "Nationality"), "ggplot")
  expect_error(hill_tuckey_pq(enterotype, "NAtionnality"))
})

data("data_fungi")
test_that("hill_tuckey_pq function works fine with data_fungi dataset", {
  expect_silent(hill_tuckey_pq(data_fungi, "Time"))
  expect_message(hill_tuckey_pq(data_fungi, "Time", silent = FALSE))
  expect_s3_class(hill_tuckey_pq(data_fungi, "Time"), "ggplot")
  expect_error(hill_tuckey_pq(data_fungi, "Timmes"))
})
