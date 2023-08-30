data(GlobalPatterns)

test_that("plot_edgeR_pq works", {
  expect_message(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), color_tax = "Kingdom"), "Perform edgeR binary test") # nolint: line_length_linter.
  expect_message(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), color_tax = "Species"), "Perform edgeR binary test")
  expect_message(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), taxolev = "Class", color_tax = "Kingdom"), "Perform edgeR binary test")
  expect_error(plot_edgeR_pq(GlobalPatterns, c("SampleType"), taxolev = "Class", color_tax = "Kingdom"), "At least one element of given pair is not a group")
  expect_error(plot_edgeR_pq(GlobalPatterns, c("SampleType", "Soil", "Feces"), color_tax = "Samples"))
})


GP <- subset_samples(GlobalPatterns, GlobalPatterns@sam_data$SampleType %in% c("Soil", "Skin"))

res <- DESeq2::DESeq(phyloseq_to_deseq2(GP, ~SampleType), test = "Wald", fitType = "local")

test_that("plot_deseq2_pq works", {
  expect_message(DESeq2::DESeq(phyloseq_to_deseq2(GP, ~SampleType), test = "Wald", fitType = "local"), "fitting model and testing")
  expect_silent(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Kingdom"))
  ## TODO expect_silent(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), color_tax = "Kingdom", verbose = TRUE))
  expect_silent(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Kingdom"))
  expect_silent(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Kingdom", taxolev = "Class"))
  expect_silent(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Class"))
  expect_silent(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), tax_table = GP@tax_table, color_tax = "Class", alpha = 0.7))
  expect_error(plot_deseq2_pq(res, c("SampleType", "Soil", "Skyp"), tax_table = GP@tax_table, color_tax = "Kingdom"))
  expect_error(plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"), color_tax = "Class"))
})

test_that("phyloseq_to_edgeR gives the good class", {
  expect_s4_class(phyloseq_to_edgeR(GlobalPatterns, "SampleType"), "DGEList")
})
