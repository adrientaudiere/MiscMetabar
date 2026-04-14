data(data_fungi)
data(enterotype)

test_that("ggaluv_pq works", {
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p <- ggaluv_pq(data_fungi, "Height")
    expect_s3_class(p, "ggplot")
  }
})

test_that("ggscatt_pq works", {
  suppressWarnings(p <- ggaluv_pq(data_fungi, wrap_factor = "Height"))
  expect_s3_class(p, "ggplot")

  expect_error(
    p <- ggaluv_pq(
      data_fungi,
      fact = "Height",
      by_sample = TRUE,
      use_ggfittext = TRUE,
      na_remove = TRUE
    )
  )

  library(ggalluvial)
  suppressWarnings(
    p <- ggaluv_pq(
      data_fungi,
      fact = "Height",
      by_sample = TRUE,
      use_ggfittext = TRUE,
      na_remove = TRUE
    )
  )
  expect_s3_class(p, "ggplot")

  suppressWarnings(
    p <- ggaluv_pq(
      data_fungi,
      fact = "Height",
      rarefy_by_sample = TRUE,
      use_geom_label = TRUE,
      rngseed = 207706,
      na_remove = TRUE
    )
  )
  expect_s3_class(p, "ggplot")
})

test_that("umap_pq works", {
  if (requireNamespace("umap", quietly = TRUE)) {
    # Regression test for issue #134: umap branch must not emit a tibble
    # .name_repair deprecation warning when converting the layout matrix.
    expect_no_warning(result <- umap_pq(data_fungi))
    expect_s3_class(result, "tbl_df")

    suppressWarnings(result <- umap_pq(data_fungi, pkg = "uwot"))
    expect_s3_class(result, "tbl_df")
  }
})
