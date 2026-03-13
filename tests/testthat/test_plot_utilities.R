data(data_fungi)
data(data_fungi_mini)
data(enterotype)


test_that("no_legend works fine", {
  p <- ggplot(mtcars) +
    geom_point(aes(x = mpg, y = hp, color = factor(cyl)))

  result <- p + no_legend()
  expect_s3_class(result, "ggplot")

  # Check that legend is removed
  built <- ggplot_build(result)
  expect_identical(result$theme$legend.position, "none")
})


test_that("ridges_pq works fine", {
  if (requireNamespace("ggridges")) {
    result <- ridges_pq(
      data_fungi_mini,
      "Time",
      alpha = 0.5,
      log10trans = FALSE
    )
    expect_s3_class(result, "ggplot")
    skip_on_cran()
    result2 <- ridges_pq(
      data_fungi_mini,
      "Time",
      alpha = 0.5,
      log10trans = TRUE
    )
    expect_s3_class(result2, "ggplot")

    result3 <- ridges_pq(data_fungi_mini, "Time", type = "ecdf")
    expect_s3_class(result3, "ggplot")
  }
})


test_that("ridges_sam_pq works fine", {
  if (requireNamespace("ggridges")) {
    result <- ridges_sam_pq(
      data_fungi_mini,
      "Time",
      alpha = 0.5,
      log10trans = FALSE
    )
    expect_s3_class(result, "ggplot")
    skip_on_cran()
    result2 <- ridges_sam_pq(
      data_fungi_mini,
      "Time",
      alpha = 0.5,
      log10trans = TRUE
    )
    expect_s3_class(result2, "ggplot")

    result3 <- ridges_sam_pq(data_fungi_mini, "Time", type = "ecdf")
    expect_s3_class(result3, "ggplot")

    result4 <- ridges_sam_pq(data_fungi_mini, "Time", nb_seq = FALSE)
    expect_s3_class(result4, "ggplot")
  }
})


test_that("treemap_pq works fine", {
  if (requireNamespace("treemapify")) {
    result2 <- suppressWarnings(treemap_pq(data_fungi_mini, "Class", "Order"))
    expect_s3_class(result2, "ggplot")
  }
})


test_that("plot_tax_pq works fine", {
  result <- plot_tax_pq(data_fungi_mini, "Height")
  expect_s3_class(result, "ggplot")
  skip_on_cran()
  result2 <- plot_tax_pq(data_fungi_mini, "Height", type = "nb_taxa")
  expect_s3_class(result2, "ggplot")

  result3 <- plot_tax_pq(data_fungi_mini, "Height", type = "both")
  expect_type(result3, "list")
  expect_length(result3, 2)
})


test_that("biplot_pq works fine", {
  data_fungi_2Height <- subset_samples(
    data_fungi_mini,
    Height %in% c("Low", "High")
  )
  skip_on_cran()
  result2 <- biplot_pq(
    data_fungi_2Height,
    fact = "Height",
    merge_sample_by = "Height"
  )
  expect_s3_class(result2, "ggplot")
})


test_that("hill_curves_pq works fine", {
  if (requireNamespace("vegan")) {
    result <- hill_curves_pq(data_fungi_mini, merge_sample_by = "Time")
    expect_s3_class(result, "ggplot")
    skip_on_cran()
    result2 <- hill_curves_pq(
      data_fungi_mini,
      color_fac = "Time",
      plot_legend = FALSE
    )
    expect_s3_class(result2, "ggplot")
  }
})


test_that("plot_edgeR_pq works fine", {
  if (requireNamespace("edgeR")) {
    data_mini_sub <- subset_samples(data_fungi_mini, !is.na(Height))
    skip_on_cran()
    result <- suppressWarnings(plot_edgeR_pq(
      data_mini_sub,
      contrast = c("Height", "Low", "High")
    ))
    expect_s3_class(result, "ggplot")
  }
})
