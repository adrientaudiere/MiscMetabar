data("data_fungi")
data("data_fungi_sp_known")
data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns
data_fungi_2trees <- subset_samples(data_fungi, data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X"))
GP_archae <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")

test_that("biplot_pq works", {
  expect_message(biplot_pq(data_fungi_2trees, merge_sample_by = "Tree_name"))
  expect_message(biplot_pq(data_fungi_2trees, fact = "Tree_name", merge_sample_by = "Tree_name"))
  expect_error(biplot_pq(data_fungi, merge_sample_by = "Tree_name"), "biplot_pq needs only two samples")
  expect_error(biplot_pq(data_fungi_2trees, fact = "Tree_name"), "biplot_pq needs only two samples")
  expect_error(biplot_pq(data_fungi_2trees, merge_sample_by = "tRREE_name"))
})

data_fungi_abun <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)
test_that("multi_biplot_pq works with data_fungi dataset", {
  p1 <- multi_biplot_pq(data_fungi_abun, split_by = "Time", na_remove = FALSE)
  p2 <- multi_biplot_pq(data_fungi_abun, "Height")
  expect_s3_class(p1[[1]], "ggplot")
  expect_type(p1, "list")
  expect_s3_class(p2[[1]], "ggplot")
  expect_type(p2, "list")
})


test_that("circle_pq works", {
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", nproc = 4))
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", min_prop_tax = 0.0001, min_prop_mod = 0.0001, nproc = 4))
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", nproc = 4, add_nb_seq = FALSE))
  expect_message(circle_pq(data_fungi_2trees, fact = "Tree_name", nproc = 4, rarefy = TRUE))
  expect_error(circle_pq(data_fungi_2trees, fact = "tRREE_name"))
})

test_that("graph_test_pq works", {
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name"))
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name", na_remove = TRUE))
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name", return_plot = FALSE))
  expect_message(graph_test_pq(subset_samples(data_fungi, !is.na(data_fungi@sam_data$Time)), fact = "Time", merge_sample_by = "Tree_name"))
  expect_error(graph_test_pq(data_fungi, fact = "Height"))
  expect_error(graph_test_pq(enterotype, fact = "Enterotype"))
  expect_error(graph_test_pq(data_fungi, fact = "tRREE_name"))
})

data_fungi <- subset_samples(data_fungi, !is.na(Time))
res_mt <- mt(data_fungi, "Time", method = "fdr", test = "f", B = 300)
test_that("plot_mt works", {
  expect_s3_class(res_mt, "data.frame")
  expect_s3_class(suppressWarnings(plot_mt(res_mt)), "ggplot")
  expect_s3_class(suppressWarnings(plot_mt(res_mt, taxa = "Genus", color_tax = "Order")), "ggplot")
})

test_that("sankey_pq works with GlobalPatterns dataset", {
  expect_silent(sankey_pq(GP))
  expect_s3_class(sankey_pq(GP), "htmlwidget")
  expect_s3_class(sankey_pq(GP), "sankeyNetwork")
  expect_silent(suppressWarnings(sankey_pq(GP, fact = "SampleType")))
  expect_silent(sankey_pq(GP, taxa = c(1:4), min_prop_tax = 0.01, units = "sequences"))
  expect_silent(sankey_pq(GP, taxa = c(1:4), min_prop_tax = 0.01, add_nb_seq = TRUE))
  expect_silent(sankey_pq(GP, taxa = c(1:4), min_prop_tax = 0.001, add_nb_seq = TRUE, tax2remove = "NRP-J"))
  expect_silent(sankey_pq(GP, taxa = c(1:4), min_prop_tax = 0.01, add_nb_seq = TRUE, units = "sequences", symbol2sub = NULL))
  expect_warning(sankey_pq(GP, taxa = c(1:4), min_prop_tax = 0.01, add_nb_seq = TRUE, units = "sequences", symbol2sub = NA))
  expect_error(sankey_pq(GP, taxa = c(1:9)))
})

test_that("sankey_pq works with data_fungi dataset", {
  expect_silent(sankey_pq(data_fungi))
  expect_s3_class(sankey_pq(data_fungi), "htmlwidget")
  expect_s3_class(sankey_pq(data_fungi), "sankeyNetwork")
  expect_silent(suppressWarnings(sankey_pq(data_fungi, fact = "Height")))
  expect_silent(sankey_pq(data_fungi, taxa = c(3:7), min_prop_tax = 0.01, units = "sequences"))
  expect_silent(sankey_pq(data_fungi, taxa = c(1:4), min_prop_tax = 0.01, add_nb_seq = TRUE))
  expect_silent(sankey_pq(data_fungi, taxa = c(1:4), min_prop_tax = 0.001, add_nb_seq = TRUE, tax2remove = "Undefined"))
  expect_silent(sankey_pq(data_fungi, taxa = c(1:4), add_nb_seq = TRUE, units = "sequences", symbol2sub = NULL))
  expect_warning(sankey_pq(data_fungi, taxa = c(1:4), min_prop_tax = 0.01, add_nb_seq = TRUE, units = "sequences", symbol2sub = NA))
  expect_error(sankey_pq(data_fungi, "HEIGHT"))
})

test_that("venn_pq works with data_fungi dataset", {
  library("grid")
  expect_silent(venn_pq(data_fungi, "Height"))
  expect_silent(venn_pq(data_fungi, "Height", min_nb_seq = 10))
  expect_silent(venn_pq(data_fungi, "Height", print_values = FALSE))
  expect_silent(venn_pq(data_fungi, "Height", print_values = FALSE) + scale_fill_hue())
  expect_silent(venn_pq(data_fungi, "Height", print_values = TRUE) + scale_fill_hue())
  expect_error(venn_pq(data_fungi))
  expect_type(venn_pq(data_fungi, "Height", print_values = TRUE) + scale_fill_hue(), "NULL")
  expect_s3_class(venn_pq(data_fungi, "Height", print_values = FALSE) + scale_fill_hue(), "ggplot")
})

test_that("ggvenn_pq works with data_fungi dataset", {
  expect_silent(ggvenn_pq(data_fungi, "Height"))
  expect_silent(suppressMessages(ggvenn_pq(data_fungi, "Height", rarefy_nb_seqs = TRUE)))
  expect_silent(ggvenn_pq(data_fungi, "Height", min_nb_seq = 2))
  expect_silent(ggvenn_pq(data_fungi, "Height", taxonomic_rank = 4))
  expect_silent(suppressMessages(ggvenn_pq(data_fungi, "Height", split_by = "Time")))
  expect_error(ggvenn_pq(data_fungi))
  expect_s3_class(ggvenn_pq(data_fungi, "Height"), "ggplot")
})

test_that("multiplot works fine", {
  res_venn1 <- ggvenn_pq(data_fungi, "Height")
  res_venn2 <- ggvenn_pq(data_fungi, "Time")
  expect_silent(multiplot(res_venn1, res_venn2))
  expect_type(multiplot(res_venn1, res_venn2), "NULL")
})

test_that("hill_pq works with data_fungi dataset", {
  expect_message(expect_message(hill_pq(data_fungi, "Height")))
  expect_message(expect_message(hill_pq(data_fungi, "Height", add_points = TRUE)))
  expect_silent(suppressMessages(hill_pq(clean_pq(subset_samples_pq(data_fungi, !is.na(data_fungi@sam_data$Height))), "Height", letters = TRUE)))
  expect_silent(suppressMessages(hill_pq(data_fungi, "Height", add_points = TRUE, color_fac = "Time")))
  expect_equal(length(hill_pq(data_fungi, "Height", add_points = TRUE)), 4)
  expect_s3_class(hill_pq(data_fungi, "Height", add_points = TRUE)[[1]], "ggplot")
})

test_that("hill_pq works with GP dataset", {
  expect_message(hill_pq(GP, "SampleType"))
  expect_message(hill_pq(GP, "SampleType", add_points = TRUE))
  expect_silent(suppressMessages(hill_pq(GP, "SampleType", letters = TRUE)))
  expect_silent(suppressMessages(hill_pq(GP, "SampleType", add_points = TRUE)))
  expect_equal(length(hill_pq(GP, "SampleType", add_points = TRUE)), 4)
  expect_s3_class(hill_pq(GP, "SampleType", add_points = TRUE)[[1]], "ggplot")
})


test_that("summary_plot_pq works with data_fungi dataset", {
  expect_message(summary_plot_pq(data_fungi))
  expect_s3_class(summary_plot_pq(data_fungi), "ggplot")
  expect_message(summary_plot_pq(data_fungi, add_info = FALSE))
  expect_message(summary_plot_pq(data_fungi, add_info = FALSE, min_seq_samples = 33))
  expect_message(summary_plot_pq(data_fungi) + scale_fill_viridis_d())
  expect_silent(summary_plot_pq(data_fungi, clean_pq = FALSE))
})

test_that("summary_plot_pq works with GP dataset", {
  expect_message(summary_plot_pq(GP))
  expect_s3_class(summary_plot_pq(GP), "ggplot")
  expect_message(summary_plot_pq(GP, add_info = FALSE))
  expect_message(summary_plot_pq(GP, add_info = FALSE, min_seq_samples = 33))
  expect_message(summary_plot_pq(GP) + scale_fill_viridis_d())
  expect_silent(summary_plot_pq(GP, clean_pq = FALSE))
})

test_that("summary_plot_pq works with enterotype dataset", {
  expect_message(summary_plot_pq(enterotype))
  expect_s3_class(summary_plot_pq(enterotype), "ggplot")
  expect_message(summary_plot_pq(enterotype, add_info = FALSE))
  expect_message(summary_plot_pq(enterotype, add_info = FALSE, min_seq_samples = 33))
  expect_message(summary_plot_pq(enterotype) + scale_fill_viridis_d())
  expect_silent(summary_plot_pq(enterotype, clean_pq = FALSE))
})

test_that("rotl_pq works with data_fungi dataset", {
  library("rotl")
  expect_s3_class(suppressWarnings(tr <- rotl_pq(data_fungi, species_colnames = "Genus_species")), "phylo")
  expect_s3_class(suppressWarnings(rotl_pq(data_fungi, species_colnames = "Genus_species", context_name = "Ascomycetes")), "phylo")
  expect_silent(plot(tr))
})

data_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
test_that("heat_tree_pq works with data_fungi dataset", {
  expect_silent(suppressMessages(ht <- heat_tree_pq(data_basidio)))
  expect_s3_class(ht, "ggplot")
})


GPsubset <- subset_taxa(
  GlobalPatterns,
  GlobalPatterns@tax_table[, 1] == "Bacteria"
)
test_that("heat_tree_pq works with GlobalPatterns dataset", {
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


test_that("accu_plot works with GlobalPatterns dataset", {
  expect_silent(accu_plot(GP_archae, fact = "X.SampleID", by.fact = TRUE))
  expect_silent(accu_plot(GP_archae, fact = "X.SampleID", by.fact = TRUE, print_sam_names = TRUE))
  expect_silent(accu_plot(GP_archae, "SampleType", add_nb_seq = TRUE, by.fact = TRUE))
  expect_silent(suppressWarnings(accu_plot(GP_archae, "SampleType", add_nb_seq = FALSE, by.fact = TRUE)))
  expect_error(accu_plot(GP_archae))
})

test_that("accu_plot works with data_fungi dataset", {
  expect_silent(accu_plot(data_fungi, fact = "Height", by.fact = TRUE))
  expect_error(accu_plot(data_fungi, fact = "Height", by.fact = FALSE))
  expect_silent(accu_plot(data_fungi, fact = "Height", by.fact = TRUE, print_sam_names = TRUE))
  expect_silent(accu_plot(data_fungi, "Height", add_nb_seq = TRUE, by.fact = TRUE))
  expect_silent(accu_plot(data_fungi, "Height", add_nb_seq = FALSE, by.fact = TRUE))
  expect_error(accu_plot(data_fungi))
})
