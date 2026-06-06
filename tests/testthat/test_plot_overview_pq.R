data(data_fungi_mini)

subset_pq <- function() {
  sn <- sample_names(data_fungi_mini)
  hi <- sn[which(data_fungi_mini@sam_data$Height == "High")[1:3]]
  lo <- sn[which(data_fungi_mini@sam_data$Height == "Low")[1:3]]
  clean_pq(prune_samples(c(hi, lo), data_fungi_mini))
}

test_that("plot_overview_pq returns a named list of ggplots by default", {
  skip_on_cran()
  skip_if_not_installed("vegan")
  skip_if_not_installed("patchwork")
  ps <- subset_pq()
  res <- suppressMessages(suppressWarnings(
    plot_overview_pq(
      ps,
      fact = "Height",
      q = c(0, 1),
      add_venn = FALSE,
      add_umap = FALSE
    )
  ))
  expect_type(res, "list")
  expect_true("ordination" %in% names(res))
  expect_true(all(vapply(res, ggplot2::is.ggplot, logical(1))))
})

test_that("plot_overview_pq assembles a patchwork when one_plot = TRUE", {
  skip_on_cran()
  skip_if_not_installed("patchwork")
  skip_if_not_installed("vegan")
  ps <- subset_pq()
  p <- suppressMessages(suppressWarnings(
    plot_overview_pq(
      ps,
      fact = "Height",
      q = c(0, 1),
      add_venn = FALSE,
      add_umap = FALSE,
      one_plot = TRUE
    )
  ))
  expect_s3_class(p, "patchwork")
})

test_that("plot_overview_pq uses a Hill scatter and skips venn for numeric fact", {
  skip_on_cran()
  skip_if_not_installed("vegan")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggstatsplot")
  ps <- subset_pq()
  expect_message(
    res <- suppressWarnings(
      plot_overview_pq(
        ps,
        fact = "Time",
        q = c(0, 1),
        add_umap = FALSE
      )
    ),
    "numeric"
  )
  expect_true("alpha" %in% names(res))
  expect_false(any(c("venn", "upset") %in% names(res)))
})

test_that("plot_overview_pq draws a Venn for a factor within venn_max", {
  skip_on_cran()
  skip_if_not_installed("vegan")
  skip_if_not_installed("ggVennDiagram")
  ps <- subset_pq()
  res <- suppressMessages(suppressWarnings(
    plot_overview_pq(
      ps,
      fact = "Height",
      q = c(0, 1),
      add_alpha = FALSE,
      add_ordination = FALSE,
      add_umap = FALSE
    )
  ))
  expect_true("venn" %in% names(res))
})

test_that("plot_overview_pq switches to UpSet above venn_max", {
  skip_on_cran()
  skip_if_not_installed("vegan")
  skip_if_not_installed("ComplexUpset")
  ps <- subset_pq()
  res <- suppressMessages(suppressWarnings(
    plot_overview_pq(
      ps,
      fact = "Height",
      q = c(0, 1),
      venn_max = 1,
      add_alpha = FALSE,
      add_ordination = FALSE,
      add_umap = FALSE
    )
  ))
  expect_true("upset" %in% names(res))
})

test_that("plot_overview_pq skips UMAP gracefully on too-few samples", {
  skip_on_cran()
  skip_if_not_installed("vegan")
  skip_if_not_installed("patchwork")
  ps <- subset_pq()
  expect_message(
    res <- suppressWarnings(
      plot_overview_pq(
        ps,
        fact = "Height",
        q = c(0, 1),
        add_venn = FALSE
      )
    ),
    "UMAP panel is skipped"
  )
  expect_false("umap" %in% names(res))
  expect_true("ordination" %in% names(res))
})

test_that("plot_overview_pq aborts on an invalid factor", {
  skip_on_cran()
  data(data_fungi_mini)
  expect_error(
    plot_overview_pq(data_fungi_mini, fact = "not_a_var"),
    "sam_data"
  )
})
