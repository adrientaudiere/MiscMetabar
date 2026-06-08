test_that("plot_volcano_pq returns a ggplot with default columns", {
  set.seed(1)
  res <- data.frame(
    log2FoldChange = rnorm(50, sd = 2),
    padj = runif(50)^3,
    Genus = sample(paste0("G", 1:10), 50, replace = TRUE)
  )
  p <- plot_volcano_pq(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano_pq handles padj == 0 and NA padj", {
  res <- data.frame(
    log2FoldChange = c(3, -3, 0.1, 2),
    padj = c(0, 1e-10, NA, 0.001)
  )
  p <- suppressMessages(plot_volcano_pq(res, lfc_threshold = 1))
  expect_s3_class(p, "ggplot")
  expect_true(all(is.finite(p$data$.neglog10)))
  # NA padj -> NotDA
  expect_identical(as.character(p$data$.status[3]), "NotDA")
  # positive, significant, large fc -> Up ; negative -> Down
  expect_identical(as.character(p$data$.status[1]), "Up")
  expect_identical(as.character(p$data$.status[2]), "Down")
})

test_that("plot_volcano_pq labels significant points when label_col is set", {
  set.seed(2)
  res <- data.frame(
    log2FoldChange = c(5, -5, 0.2),
    padj = c(1e-5, 1e-4, 0.9),
    Genus = c("A", "B", "C")
  )
  p <- plot_volcano_pq(res, label_col = "Genus", label_n = 2)
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano_pq accepts custom fc/padj column names", {
  res <- data.frame(effect = c(2, -2), wi.eBH = c(0.01, 0.02))
  p <- plot_volcano_pq(res, fc = "effect", padj = "wi.eBH")
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano_pq aborts on missing columns", {
  res <- data.frame(a = 1:3, b = 1:3)
  expect_error(plot_volcano_pq(res), "not found")
})

test_that("plot_volcano_pq auto-detects ALDEx2 columns", {
  res <- data.frame(
    effect = c(2, -2, 0.1),
    wi.eBH = c(0.01, 0.03, 0.8),
    we.eBH = c(0.02, 0.04, 0.9)
  )
  p <- plot_volcano_pq(res)
  expect_s3_class(p, "ggplot")
  expect_identical(as.character(p$data$.status[1]), "Up")
  expect_identical(as.character(p$data$.status[2]), "Down")
})

test_that("plot_volcano_pq auto-detects ANCOMBC list and extracts $res", {
  ancombc_res <- list(
    res = data.frame(
      taxon = letters[1:4],
      lfc_GroupB = c(2, -2, 0.1, -0.1),
      q_GroupB = c(0.01, 0.02, 0.8, 0.9)
    )
  )
  p <- suppressMessages(plot_volcano_pq(ancombc_res))
  expect_s3_class(p, "ggplot")
  expect_identical(as.character(p$data$.status[1]), "Up")
  expect_identical(as.character(p$data$.status[2]), "Down")
})

test_that("plot_volcano_pq ANCOMBC: explicit fc/padj overrides auto-detect", {
  ancombc_res <- list(
    res = data.frame(
      taxon = letters[1:3],
      lfc_GroupB = c(2, -2, 0.1),
      q_GroupB = c(0.01, 0.02, 0.8),
      lfc_GroupC = c(-1, 1, 0.3),
      q_GroupC = c(0.03, 0.04, 0.7)
    )
  )
  p <- suppressMessages(
    plot_volcano_pq(ancombc_res, fc = "lfc_GroupC", padj = "q_GroupC")
  )
  expect_s3_class(p, "ggplot")
  expect_identical(as.character(p$data$.status[2]), "Up")
})
