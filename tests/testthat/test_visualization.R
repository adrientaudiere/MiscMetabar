library("divent")

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

test_that("umap_pq caps n_neighbors on small phyloseq objects", {
  data_8 <- prune_samples(sample_names(data_fungi_mini)[1:8], data_fungi_mini)
  if (requireNamespace("umap", quietly = TRUE)) {
    # default n_neighbors (15) and a too large value are capped to 7
    expect_equal(nrow(umap_pq(data_8)), 8)
    expect_equal(nrow(umap_pq(data_8, n_neighbors = 30)), 8)
    expect_equal(nrow(umap_pq(data_8, n_neighbors = 3)), 8)
  }
  if (requireNamespace("uwot", quietly = TRUE)) {
    suppressWarnings(result <- umap_pq(data_8, pkg = "uwot"))
    expect_equal(nrow(result), 8)
  }
})
