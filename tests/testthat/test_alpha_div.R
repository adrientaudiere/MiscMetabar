data(data_fungi)
data(data_fungi_mini)


test_that("hill_tuckey_pq works with different parameters", {
  data("GlobalPatterns", package = "phyloseq")
  GlobalPatterns@sam_data[, "Soil_logical"] <-
    ifelse(
      GlobalPatterns@sam_data[, "SampleType"] == "Soil",
      "Soil",
      "Not Soil"
    )

  expect_s3_class(
    suppressMessages(hill_tuckey_pq(GlobalPatterns, "Soil_logical")),
    "ggplot"
  )
  skip_on_cran()
  expect_s3_class(
    suppressMessages(hill_tuckey_pq(
      GlobalPatterns,
      "Soil_logical",
      hill_scales = 1:2
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressMessages(hill_tuckey_pq(
      GlobalPatterns,
      "Soil_logical",
      correction_for_sample_size = FALSE
    )),
    "ggplot"
  )
})


test_that("hill_test_rarperm_pq works with data_fungi", {
  if (requireNamespace("ggstatsplot")) {
    skip_on_cran()
    expect_type(
      suppressWarnings(suppressMessages(hill_test_rarperm_pq(
        data_fungi,
        "Time",
        nperm = 2
      ))),
      "list"
    )
    result <- suppressWarnings(suppressMessages(
      hill_test_rarperm_pq(data_fungi, "Height", nperm = 3, p_val_signif = 0.9)
    ))
    expect_type(result, "list")
    expect_true("method" %in% names(result))
    expect_true("expressions" %in% names(result))
    expect_true("plots" %in% names(result))
    expect_true("pvals" %in% names(result))
    expect_true("prop_signif" %in% names(result))
    expect_true("statistics" %in% names(result))
  }
})


test_that("glmutli_pq works with data_fungi", {
  skip_on_cran()
  if (requireNamespace("glmulti")) {
    result <- suppressWarnings(glmutli_pq(
      data_fungi_mini,
      "Hill_0 ~ Abundance + Time",
      level = 1
    ))
    expect_s3_class(result, "data.frame")
    expect_true("estimates" %in% colnames(result))
    expect_true("importance" %in% colnames(result))
  }
})


test_that("glmutli_pq works with more complex scheme", {
  if (requireNamespace("glmulti", quietly = TRUE)) {
    res_glmulti <-
      glmutli_pq(
        data_fungi,
        "Hill_0 ~ Hill_1 + Abundance + Time + Height",
        level = 1
      )
    expect_equal(dim(res_glmulti), c(5, 6))
    res_glmulti_interaction <-
      glmutli_pq(data_fungi, "Hill_0 ~ Abundance + Time + Height", level = 2)
    expect_equal(dim(res_glmulti_interaction), c(11, 6))
  }
})
