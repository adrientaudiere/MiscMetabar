data(data_fungi)
data(enterotype)

test_that("ggaluv_pq works", {
  skip_on_cran()
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p <- ggaluv_pq(data_fungi, "Height")
    expect_s3_class(p, "ggplot")
  }
})

test_that("ggscatt_pq works", {
  skip_on_cran()
  p <- ggscatt_pq(data_fungi, "Height")
  expect_s3_class(p, "ggplot")
})

test_that("umap_pq works", {
  skip_on_cran()
  if (requireNamespace("umap", quietly = TRUE)) {
    result <- umap_pq(data_fungi)
    expect_s4_class(result, "phyloseq")
  }
})
