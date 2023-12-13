data("GlobalPatterns", package = "phyloseq")
GP <- subset_samples(GlobalPatterns, GlobalPatterns@sam_data$SampleType %in% c("Soil", "Skin"))

fac2col <- function(x, col.pal = funky_color, na.col = "transparent", seed = NULL) {
  x <- factor(x)
  lev <- levels(x)
  nlev <- length(lev)
  if (!is.null(seed)) {
    set.seed(seed)
    newseed <- round(runif(1, 1, 1e+09))
    on.exit(set.seed(newseed))
    col <- sample(col.pal(nlev))
  } else {
    col <- col.pal(nlev)
  }
  res <- rep(na.col, length(x))
  res[!is.na(x)] <- col[as.integer(x[!is.na(x)])]
  return(res)
}

test_that("plot_edgeR_pq works", {
  expect_message(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), color_tax = "Kingdom"), "Perform edgeR binary test") # nolint: line_length_linter.
  expect_message(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), color_tax = "Species"), "Perform edgeR binary test")
  expect_message(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), taxolev = "Class", color_tax = "Kingdom"), "Perform edgeR binary test")
  expect_error(plot_edgeR_pq(GlobalPatterns, "SampleType", taxolev = "Class", color_tax = "Kingdom"), "At least one element of given pair is not a group")
  expect_error(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), color_tax = "Samples"))
})
test_that("plot_deseq2_pq works with results", {
  expect_message(res <- DESeq2::DESeq(phyloseq_to_deseq2(GP, ~SampleType), test = "Wald", fitType = "local"), "fitting model and testing")
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Kingdom")))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, tax_depth = "Genus")))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"),
    tax_table = GP@tax_table, tax_depth = "Family",
    color_tax = fac2col(as.vector(GP@tax_table[, "Order"]))
  )))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Kingdom", verbose = TRUE)))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Kingdom")))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Kingdom", taxolev = "Class")))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Class")))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Class", alpha = 0.7)))
  expect_error(plot_deseq2_pq(res, c("SampleType", "Soil", "Skyp"), tax_table = GP@tax_table, color_tax = "Kingdom"))
  expect_error(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), color_tax = "Class"))
  expect_message(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Class", select_taxa = "522457"))
  expect_message(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Class", select_taxa = c("522457", "271582")))
  expect_silent(suppressWarnings(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, select_taxa = c("522457", "200359"))))
})

test_that("plot_deseq2_pq works with phyloseq object", {
  expect_message(suppressWarnings(plot_deseq2_pq(GP, c("SampleType", "Soil", "Skin"))))
  expect_message(suppressWarnings(plot_deseq2_pq(GP, c("SampleType", "Soil", "Skin"), color_tax = "Class", select_taxa = c("522457", "271582", "200359"))))
  expect_message(suppressWarnings(plot_deseq2_pq(GP, c("SampleType", "Soil", "Skin"), taxolev = "Class", verbose = TRUE)))
  expect_message(suppressWarnings(plot_deseq2_pq(GP, c("SampleType", "Soil", "Skin"), alpha = 0.1)))
})

GlobalPatterns_row <- clean_pq(GP, force_taxa_as_columns = TRUE)

test_that("phyloseq_to_edgeR gives the good class", {
  expect_s4_class(phyloseq_to_edgeR(GlobalPatterns, "SampleType"), "DGEList")
  expect_s4_class(phyloseq_to_edgeR(GlobalPatterns_row, "SampleType"), "DGEList")
})
