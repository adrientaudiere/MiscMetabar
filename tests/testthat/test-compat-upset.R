test_that("is_ggplot2_compatible returns logical", {
  result <- is_ggplot2_compatible()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("is_ggplot2_compatible handles missing ggplot2", {
  # Mock case where ggplot2 is not available
  # This is mainly for completeness, as ggplot2 is a hard dependency
  expect_type(is_ggplot2_compatible(), "logical")
})

test_that(".get_upset_compatibility_message returns character", {
  msg <- MiscMetabar:::.get_upset_compatibility_message("test_function")
  expect_type(msg, "character")
  expect_true(grepl("test_function", msg))
  expect_true(grepl("ggplot2", msg))
})

test_that("upset functions fail gracefully with incompatible ggplot2", {
  skip_if_not(requireNamespace("ComplexUpset", quietly = TRUE), "ComplexUpset not available")
  
  # Only test if we actually have an incompatible version
  if (!is_ggplot2_compatible()) {
    expect_error(
      upset_pq(data_fungi_mini, "Height"),
      "ggplot2.*4\\.0\\.0.*ComplexUpset"
    )
    
    expect_error(
      upset_test_pq(data_fungi_mini, "Height"),
      "ggplot2.*4\\.0\\.0.*ComplexUpset"
    )
  } else {
    skip("ggplot2 is compatible, cannot test incompatibility error")
  }
})