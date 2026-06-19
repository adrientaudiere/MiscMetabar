data("data_fungi_mini", package = "MiscMetabar")
df_env <- subset_samples(data_fungi_mini, !is.na(Height) & !is.na(Time))

test_that("plot_ordination_pq returns a ggplot without envfit", {
  skip_if_not_installed("vegan")
  p <- suppressWarnings(suppressMessages(plot_ordination_pq(
    df_env,
    method = "bray"
  )))
  expect_s3_class(p, "ggplot")
})

test_that("add_envfit overlays vectors for a continuous variable", {
  skip_if_not_installed("vegan")
  p0 <- suppressWarnings(suppressMessages(plot_ordination_pq(
    df_env,
    method = "bray"
  )))
  p1 <- suppressWarnings(suppressMessages(plot_ordination_pq(
    df_env,
    method = "bray",
    add_envfit = TRUE,
    envfit_fact = "Time"
  )))
  expect_s3_class(p1, "ggplot")
  expect_gt(length(p1$layers), length(p0$layers))
})

test_that("add_envfit overlays centroids for a factor variable", {
  skip_if_not_installed("vegan")
  p1 <- suppressWarnings(suppressMessages(plot_ordination_pq(
    df_env,
    method = "bray",
    add_envfit = TRUE,
    envfit_fact = "Height"
  )))
  expect_s3_class(p1, "ggplot")
})
