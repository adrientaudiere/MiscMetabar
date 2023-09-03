data("data_fungi")
data("data_fungi_sp_known")
data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns
data_fungi_2trees <- subset_samples(data_fungi, data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X"))
GP_archae <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")

test_that("rotl_pq works with data_fungi dataset", {
  library("rotl")
  expect_s3_class(suppressWarnings(tr <- rotl_pq(data_fungi, species_colnames = "Genus_species")), "phylo")
  expect_s3_class(suppressWarnings(rotl_pq(data_fungi, species_colnames = "Genus_species", context_name = "Ascomycetes")), "phylo")
  expect_silent(plot(tr))
})

data_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
test_that("heat_tree_pq works with data_fungi dataset", {
  library(metacoder)
  expect_silent(suppressMessages(ht <- heat_tree_pq(data_basidio)))
  expect_s3_class(ht, "ggplot")
})

GPsubset <- subset_taxa(
  GlobalPatterns,
  GlobalPatterns@tax_table[, 1] == "Bacteria"
)
test_that("heat_tree_pq works with GlobalPatterns dataset", {
  library(metacoder)
  expect_silent(suppressMessages(ht <- heat_tree_pq(GPsubset)))
  expect_silent(suppressMessages(ht <- heat_tree_pq(GPsubset, node_size = n_obs, node_color = n_obs, node_label = taxon_names, tree_label = taxon_names, node_size_trans = "log10 area")))
  expect_s3_class(ht, "ggplot")
})


test_that("plot_tax_pq works with data_fungi dataset", {
  expect_silent(suppressMessages(pt <- plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Class", add_info = FALSE)))
  expect_silent(suppressMessages(pt <- plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Class")))
  expect_s3_class(pt, "ggplot")
  expect_silent(suppressMessages(pt <- plot_tax_pq(data_fungi_sp_known, "Time", taxa_fill = "Class")))
  expect_s3_class(pt, "ggplot")
  expect_silent(suppressMessages(pt <- plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Class", type = "nb_asv", add_info = FALSE)))
  expect_silent(suppressMessages(pt <- plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Class", type = "nb_asv")))
  expect_s3_class(pt, "ggplot")
  expect_silent(suppressMessages(pt <- plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Class", type = "both", add_info = FALSE)))
  expect_silent(suppressMessages(pt <- plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Class", type = "both")))
  expect_s3_class(pt[[1]], "ggplot")
  expect_silent(suppressMessages(plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Class", na_remove = TRUE)))
  expect_silent(suppressMessages(plot_tax_pq(data_fungi_sp_known, "Time", merge_sample_by = "Time", taxa_fill = "Order", clean_pq = FALSE)))
  expect_silent(suppressMessages(plot_tax_pq(data_fungi_sp_known, "Height", merge_sample_by = "Height", taxa_fill = "Order")))
})

