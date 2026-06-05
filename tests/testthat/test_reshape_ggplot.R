test_that("reshape_ggplot wraps long titles and returns a ggplot", {
  skip_on_cran()
  skip_if_not_installed("stringr")
  long_title <- paste(rep("word", 30), collapse = " ")
  p <- ggplot2::ggplot(
    data.frame(x = 1:3, y = 1:3),
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(title = long_title, subtitle = long_title)
  res <- reshape_ggplot(p, width = 20)
  expect_s3_class(res, "ggplot")
  expect_match(res$labels$title, "\n")
  expect_match(res$labels$subtitle, "\n")
})

test_that("reshape_ggplot does not error when axis labels are missing", {
  skip_on_cran()
  skip_if_not_installed("stringr")
  p <- ggplot2::ggplot(
    data.frame(x = 1:3, y = 1:3),
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(x = NULL, y = NULL, title = NULL)
  expect_s3_class(reshape_ggplot(p), "ggplot")
})
