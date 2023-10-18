data(data_fungi)
data("enterotype")

test_that("dist_bycol works fine", {
  expect_equal(length(
    dist_bycol(
      data_fungi@otu_table,
      as_binary_otu_table(data_fungi)@otu_table
    )
  ), 2)
  expect_error(length(dist_bycol(
    data_fungi@otu_table, enterotype@otu_table
  )))
})

test_that("all_object_size works fine", {
  expect_type(all_object_size(), "double")
})

test_that("diff_fct_diff_class works fine", {
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Sample_id,
      numeric_fonction = sum,
      na.rm = TRUE
    ),
    17852
  )
  expect_equal(
    round(diff_fct_diff_class(
      data_fungi@sam_data$Time,
      numeric_fonction = mean,
      na.rm = TRUE
    ),2),
    5.80
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height == "Low",
      logical_method = "TRUE_if_one"
    ),
    TRUE
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height == "Low",
      logical_method = "NA_if_not_all_TRUE"
    ),
    NA
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height == "Low",
      logical_method = "FALSE_if_not_all_TRUE"
    ),
    FALSE
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height,
      character_method = "unique_or_na"
    ),
    NA
  )
  expect_equal(
    diff_fct_diff_class(
      c("IE", "IE"),
      character_method = "unique_or_na"
    ),
    "IE"
  )
  expect_equal(
    diff_fct_diff_class(
      c("IE", "IE", "TE","TE"),
      character_method = "more_frequent"
    ),
    "IE"
  )
  expect_equal(
    diff_fct_diff_class(
      c("IE", "IE", "TE","TE"),
      character_method = "more_frequent_without_equality"
    ),
    NA
  )
})
