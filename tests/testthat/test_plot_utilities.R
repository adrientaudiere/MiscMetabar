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
  expect_equal(result$theme$legend.position, "none")
})


test_that("bar_pq works fine", {
  result <- bar_pq(data_fungi_mini, fact = "Height", taxa = "Class")
  expect_s3_class(result, "ggplot")
  skip_on_cran()
  result2 <- bar_pq(data_fungi_mini, fact = "Height", taxa = "Class", percent_bar = TRUE)
  expect_s3_class(result2, "ggplot")
  
  result3 <- bar_pq(data_fungi_mini, fact = "Height", taxa = "Class", nb_seq = FALSE)
  expect_s3_class(result3, "ggplot")
})


test_that("ridges_pq works fine", {
  if (requireNamespace("ggridges")) {
    result <- ridges_pq(data_fungi_mini, "Time", alpha = 0.5, log10trans = FALSE)
    expect_s3_class(result, "ggplot")
    skip_on_cran()
    result2 <- ridges_pq(data_fungi_mini, "Time", alpha = 0.5, log10trans = TRUE)
    expect_s3_class(result2, "ggplot")
    
    result3 <- ridges_pq(data_fungi_mini, "Time", type = "ecdf")
    expect_s3_class(result3, "ggplot")
  }
})


test_that("treemap_pq works fine", {
  if (requireNamespace("treemapify")) {
    result <- suppressWarnings(treemap_pq(data_fungi_mini, "Class"))
    expect_s3_class(result, "ggplot")
    skip_on_cran()
    result2 <- suppressWarnings(treemap_pq(data_fungi_mini, "Class", subgroup = "Order"))
    expect_s3_class(result2, "ggplot")
  }
})


test_that("plot_tax_pq works fine", {
  result <- plot_tax_pq(data_fungi_mini, "Height")
  expect_type(result, "list")
  expect_s3_class(result[[1]], "ggplot")
  skip_on_cran()
  result2 <- plot_tax_pq(data_fungi_mini, "Height", type = "nb_taxa")
  expect_s3_class(result2, "ggplot")
  
  result3 <- plot_tax_pq(data_fungi_mini, "Height", type = "both")
  expect_type(result3, "list")
  expect_length(result3, 2)
})


test_that("biplot_pq works fine", {
  result <- biplot_pq(data_fungi_mini, fact = "Height")
  expect_s3_class(result, "ggplot")
  skip_on_cran()
  result2 <- biplot_pq(data_fungi_mini, fact = "Height", merge_sample_by = "Height")
  expect_s3_class(result2, "ggplot")
})


test_that("hill_curves_pq works fine", {
  if (requireNamespace("vegan")) {
    result <- hill_curves_pq(data_fungi_mini, merge_sample_by = "Time")
    expect_s3_class(result, "ggplot")
    skip_on_cran()
    result2 <- hill_curves_pq(data_fungi_mini, color_fac = "Time", plot_legend = FALSE)
    expect_s3_class(result2, "ggplot")
  }
})


test_that("plot_edgeR_pq works fine", {
  if (requireNamespace("edgeR")) {
    data_mini_sub <- subset_samples(data_fungi_mini, !is.na(Height))
    skip_on_cran()
    result <- suppressWarnings(plot_edgeR_pq(data_mini_sub,
      fact = "Height",
      levels_fact = c("Low", "High")
    ))
    expect_type(result, "list")
  }
})


test_that("ggScree works fine", {
  skip_on_cran()
  ord <- ordinate(data_fungi_mini, "PCoA", "bray")
  result <- ggScree(ord)
  expect_s3_class(result, "ggplot")
})


test_that("plot_dist_as_heatmap works fine", {
  skip_on_cran()
  dist_mat <- vegan::vegdist(t(data_fungi_mini@otu_table), method = "bray")
  result <- plot_dist_as_heatmap(dist_mat)
  expect_s3_class(result, "ggplot")
})
