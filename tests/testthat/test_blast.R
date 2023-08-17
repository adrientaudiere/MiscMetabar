data("GlobalPatterns")
GlobalPatterns@sam_data[, "Soil_logical"] <- ifelse(GlobalPatterns@sam_data[, "SampleType"] == "Soil", "Soil", "Not Soil")

test_that("write_pq function works fine with enterotype dataset", {
  expect_message(hill_tuckey_pq(GlobalPatterns))
})
