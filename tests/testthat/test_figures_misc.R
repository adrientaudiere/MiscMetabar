data("data_fungi")
data("data_fungi_sp_known")
data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns
data_fungi_2trees <- subset_samples(data_fungi, data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X"))
GP_archae <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
data_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")

test_that("tsne_pq works with data_fungi dataset", {
  expect_silent(suppressMessages(res_tsne <- tsne_pq(data_fungi)))
  expect_s3_class(res_tsne, "Rtsne")
  expect_silent(suppressMessages(res_tsne <- tsne_pq(data_fungi, dims = 3, perplexity = 25)))
})

test_that("plot_tsne_pq works with data_fungi dataset", {
  expect_silent(suppressMessages(pt <- plot_tsne_pq(data_fungi, fact = "Height", perplexity = 15)))
  expect_s3_class(pt, "ggplot")
  expect_message(plot_tsne_pq(data_fungi))
})

test_that("SRS_curve_pq works with data_fungi dataset", {
  expect_silent(suppressMessages(sc <- SRS_curve_pq(data_basidio)))
  expect_silent(suppressMessages(sc <- SRS_curve_pq(data_basidio, clean_pq = TRUE)))
  expect_s3_class(sc, "recordedplot")
  expect_silent(suppressMessages(sc <- SRS_curve_pq(data_basidio, metric = "shannon")))
  expect_silent(suppressMessages(sc <- SRS_curve_pq(data_basidio, step = 20, rarefy.repeats = 15)))
})

test_that("multiplot works fine", {
  res_venn1 <- ggvenn_pq(data_fungi, "Height")
  res_venn2 <- ggvenn_pq(data_fungi, "Time")
  expect_silent(multiplot(res_venn1, res_venn2))
  expect_message(multiplot(res_venn1))
  expect_type(multiplot(res_venn1, res_venn2), "NULL")
})
