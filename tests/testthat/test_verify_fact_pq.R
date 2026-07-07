data(data_fungi_mini)

test_that("verify_fact_pq passes for existing columns", {
  expect_invisible(verify_fact_pq(data_fungi_mini, fact = "Height"))
  expect_identical(
    verify_fact_pq(data_fungi_mini, fact = "Height"),
    data_fungi_mini
  )
  expect_no_error(verify_fact_pq(data_fungi_mini, modality = "Time"))
  expect_no_error(verify_fact_pq(data_fungi_mini, fact = c("Height", "Time")))
})

test_that("verify_fact_pq errors on a missing column", {
  expect_error(
    verify_fact_pq(data_fungi_mini, fact = "Heigth"),
    "not found"
  )
  expect_error(
    verify_fact_pq(data_fungi_mini, modality = "does_not_exist"),
    "not found"
  )
})

test_that("verify_fact_pq enforces exactly two levels for bifactor", {
  expect_error(
    verify_fact_pq(data_fungi_mini, bifactor = "Height"),
    "exactly two levels"
  )
  data_2h <- subset_samples_pq(
    data_fungi_mini,
    data_fungi_mini@sam_data$Height %in% c("Low", "High")
  )
  expect_no_error(verify_fact_pq(data_2h, bifactor = "Height"))
})
