data("data_fungi")
data("data_fungi_sp_known")
data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns
data_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
data_basidio_2trees <-
  subset_samples(data_basidio, Tree_name %in% c("A10-005", "AD30-abm-X"))
GP_archae <-
  subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")

data_fungi_woNA4time <- subset_samples(data_fungi,!is.na(Time))
res_mt <-
  mt(
    data_fungi_woNA4time,
    "Time",
    method = "fdr",
    test = "f",
    B = 300
  )

test_that("circle_pq works", {
  expect_message(circle_pq(data_basidio_2trees, fact = "Tree_name", nproc = 1))
  expect_message(circle_pq(
    clean_pq(data_basidio_2trees, force_taxa_as_rows = TRUE),
    fact = "Tree_name",
    nproc = 1
  ))
  expect_message(circle_pq(
    data_basidio,
    fact = "Tree_name",
    nproc = 1,
    min_prop_mod = 0.001
  ))
  expect_message(circle_pq(
    data_basidio_2trees,
    fact = "Tree_name",
    log10trans = TRUE,
    nproc = 1
  ))
  expect_message(
    circle_pq(
      data_basidio_2trees,
      fact = "Tree_name",
      min_prop_tax = 0.0001,
      min_prop_mod = 0.0001,
      nproc = 1
    )
  )
  expect_message(circle_pq(
    data_basidio_2trees,
    fact = "Tree_name",
    nproc = 1,
    add_nb_seq = FALSE
  ))
  expect_message(circle_pq(
    data_basidio_2trees,
    fact = "Tree_name",
    nproc = 1,
    rarefy = TRUE
  ))
  expect_error(circle_pq(data_basidio_2trees, fact = "tRREE_name"))
  expect_error(circle_pq(
    data_basidio_2trees@otu_table,
    fact = "Tree_name",
    nproc = 1
  ))
})

test_that("graph_test_pq works", {
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name"))
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name", na_remove = TRUE))
  expect_silent(graph_test_pq(data_fungi, fact = "Tree_name", return_plot = FALSE))
  expect_message(graph_test_pq(
    subset_samples(data_fungi,!is.na(data_fungi@sam_data$Time)),
    fact = "Time",
    merge_sample_by = "Tree_name"
  ))
  expect_error(graph_test_pq(data_fungi, fact = "Height"))
  expect_error(graph_test_pq(enterotype, fact = "Enterotype"))
  expect_error(graph_test_pq(data_fungi, fact = "tRREE_name"))
})


test_that("plot_mt works", {
  expect_s3_class(res_mt, "data.frame")
  expect_s3_class(suppressWarnings(plot_mt(res_mt)), "ggplot")
  expect_s3_class(suppressWarnings(plot_mt(
    res_mt,
    taxa = "Genus", color_tax = "Order"
  )), "ggplot")
})

test_that("sankey_pq works with GlobalPatterns dataset", {
  expect_silent(sankey_pq(GP_archae))
  expect_s3_class(sankey_pq(GP_archae), "htmlwidget")
  expect_s3_class(sankey_pq(GP_archae), "sankeyNetwork")
  expect_silent(suppressWarnings(sankey_pq(GP_archae, fact = "SampleType")))
  expect_silent(sankey_pq(
    GP_archae,
    taxa = 1:4,
    min_prop_tax = 0.01,
    units = "sequences"
  ))
  expect_silent(sankey_pq(
    GP_archae,
    taxa = 1:4,
    min_prop_tax = 0.01,
    add_nb_seq = TRUE
  ))
  expect_silent(suppressWarnings(sankey_pq(
    GP_archae,
    fact = "SampleType", add_nb_seq = TRUE
  )))
  expect_silent(
    sankey_pq(
      GP_archae,
      taxa = 1:4,
      min_prop_tax = 0.001,
      add_nb_seq = TRUE,
      tax2remove = "NRP-J"
    )
  )
  expect_silent(
    sankey_pq(
      GP_archae,
      taxa = 1:4,
      min_prop_tax = 0.01,
      add_nb_seq = TRUE,
      units = "sequences",
      symbol2sub = NULL
    )
  )
  expect_warning(
    sankey_pq(
      GP_archae,
      taxa = 1:4,
      min_prop_tax = 0.01,
      add_nb_seq = TRUE,
      units = "sequences",
      symbol2sub = NA
    )
  )
  expect_error(sankey_pq(GP_archae, taxa = 1:9))
  expect_error(sankey_pq(GP_archae@otu_table))
})

test_that("sankey_pq works with data_fungi dataset", {
  expect_silent(sankey_pq(data_fungi))
  expect_s3_class(sankey_pq(data_fungi), "htmlwidget")
  expect_s3_class(sankey_pq(data_fungi), "sankeyNetwork")
  expect_silent(suppressWarnings(sankey_pq(data_fungi, fact = "Height")))
  expect_silent(sankey_pq(
    data_fungi,
    taxa = 3:7,
    min_prop_tax = 0.01,
    units = "sequences"
  ))
  expect_silent(sankey_pq(
    data_fungi,
    taxa = 1:4,
    min_prop_tax = 0.01,
    add_nb_seq = TRUE
  ))
  expect_silent(
    sankey_pq(
      data_fungi,
      taxa = 1:4,
      min_prop_tax = 0.001,
      add_nb_seq = TRUE,
      tax2remove = "Undefined"
    )
  )
  expect_silent(
    sankey_pq(
      data_fungi,
      taxa = 1:4,
      add_nb_seq = TRUE,
      units = "sequences",
      symbol2sub = NULL
    )
  )
  expect_warning(
    sankey_pq(
      data_fungi,
      taxa = 1:4,
      min_prop_tax = 0.01,
      add_nb_seq = TRUE,
      units = "sequences",
      symbol2sub = NA
    )
  )
  expect_error(sankey_pq(data_fungi, "HEIGHT"))
})

test_that("venn_pq works with data_fungi dataset", {
  library("grid")
  expect_silent(venn_pq(data_fungi, "Height"))
  expect_silent(suppressMessages(venn_pq(
    clean_pq(
      data_fungi,
      force_taxa_as_rows = TRUE,
      silent = TRUE
    ),
    "Height"
  )))
  expect_silent(venn_pq(data_fungi, "Height", min_nb_seq = 10))
  expect_silent(venn_pq(data_fungi, "Height", print_values = FALSE))
  expect_silent(venn_pq(data_fungi, "Height", print_values = FALSE) + scale_fill_hue())
  expect_silent(venn_pq(data_fungi, "Height", print_values = TRUE) + scale_fill_hue())
  expect_error(venn_pq(data_fungi))
  expect_error(venn_pq(data_fungi@otu_table, "Height"))
  expect_type(venn_pq(data_fungi, "Height", print_values = TRUE) + scale_fill_hue(),
              "NULL")
  expect_s3_class(venn_pq(data_fungi, "Height", print_values = FALSE) + scale_fill_hue(),
                  "ggplot")
})

test_that("ggvenn_pq works with data_fungi dataset", {
  expect_silent(ggvenn_pq(data_fungi, "Height"))
  expect_message(ggvenn_pq(data_fungi, "Height", rarefy_before_merging = TRUE))
  expect_message(ggvenn_pq(data_fungi, "Height", rarefy_after_merging = TRUE))
  expect_silent(ggvenn_pq(data_fungi, "Height", add_nb_sequences = TRUE))
  expect_silent(suppressMessages(ggvenn_pq(data_fungi, "Height", rarefy_nb_seqs = TRUE)))
  expect_silent(ggvenn_pq(data_fungi, "Height", min_nb_seq = 2))
  expect_silent(ggvenn_pq(data_fungi, "Height", taxonomic_rank = 4))
  expect_silent(suppressMessages(ggvenn_pq(data_fungi, "Height", split_by = "Time")))
  expect_error(ggvenn_pq(data_fungi))
  expect_s3_class(ggvenn_pq(data_fungi, "Height"), "ggplot")
  expect_error(ggvenn_pq(data_fungi@otu_table, "Height"))
})


test_that("upset_pq works with data_fungi dataset", {
  expect_silent(suppressMessages(upset_pq(data_fungi, "Height")))
  expect_s3_class(upset_pq(data_fungi, "Height", taxa_fill = "Class"), "ggplot")
  expect_s3_class(upset_pq(data_fungi, "Height"), "ggplot")
  expect_s3_class(
    upset_pq(
      data_fungi,
      "Height",
      na_remove = TRUE,
      rarefy_after_merging = TRUE
    ),
    "ggplot"
  )

  expect_s3_class(upset_pq(data_fungi, "Time"), "ggplot")
  expect_s3_class(upset_pq(data_fungi, "Time", min_nb_seq = 10), "ggplot")
  expect_s3_class(upset_pq(data_fungi, "Time", numeric_fonction = mean),
                  "ggplot")

  expect_error(upset_pq(data_fungi))
})

test_that("upset_test_pq works with data_fungi dataset", {
  expect_s3_class(upset_test_pq(data_fungi, "Height"), "data.frame")
  expect_s3_class(upset_test_pq(data_fungi, "Time"), "data.frame")
  expect_s3_class(upset_test_pq(data_fungi, "Time", min_nb_seq = 10),
                  "data.frame")
  expect_s3_class(upset_test_pq(data_fungi, "Time", numeric_fonction = mean),
                  "data.frame")
  expect_s3_class(
    upset_test_pq(
      data_fungi,
      "Time",
      numeric_fonction = mean,
      var_to_test = c("OTU", "Guild", "Genus")
    ),
    "data.frame"
  )
  expect_error(upset_test_pq(data_fungi, "Height", var_to_test = "GUILDDDS"))
  expect_error(upset_test_pq(data_fungi))
})

test_that("plot_LCBD_pq works with data_fungi dataset", {
  expect_s3_class(plot_LCBD_pq(
    data_fungi,
    nperm = 100,
    only_plot_significant = FALSE
  ),
  "ggplot")
  expect_s3_class(
    plot_LCBD_pq(
      data_fungi,
      nperm = 100,
      only_plot_significant = TRUE,
      pval = 0.2
    ),
    "ggplot"
  )
  expect_s3_class(
    plot_LCBD_pq(
      data_fungi,
      nperm = 100,
      only_plot_significant = TRUE,
      p_adjust_method = "holm",
      sam_variables = c("Time", "Height")
    ),
    "ggplot"
  )
})


test_that("LCBD_pq works with data_fungi dataset", {
  expect_s3_class(LCBD_pq(data_fungi, nperm = 100), "beta.div")
  expect_s3_class(LCBD_pq(data_fungi, nperm = 100, method = "jaccard"),
                  "beta.div")
})

test_that("plot_LCBD_pq works with data_fungi dataset", {
  expect_s3_class(plot_LCBD_pq(
    data_fungi,
    nperm = 100,
    only_plot_significant = FALSE
  ),
  "ggplot")
  expect_s3_class(
    plot_LCBD_pq(
      data_fungi,
      nperm = 100,
      only_plot_significant = TRUE,
      pval = 0.2
    ),
    "ggplot"
  )
  expect_s3_class(
    plot_LCBD_pq(
      data_fungi,
      nperm = 100,
      only_plot_significant = TRUE,
      p_adjust_method = "holm",
      sam_variables = c("Time", "Height")
    ),
    "ggplot"
  )
})


test_that("plot_SCBD_pq works with data_fungi dataset", {
  expect_s3_class(plot_SCBD_pq(data_fungi), "ggplot")
  expect_s3_class(
    plot_SCBD_pq(
      data_fungi,
      tax_level = "Class",
      tax_color = "Phylum",
      min_SCBD = 0
    ) ,
    "ggplot"
  )
})

test_that("multipatt_pq works with data_fungi dataset", {
  expect_s3_class(multipatt_pq(subset_samples(data_fungi, !is.na(Time)),
                               fact = "Time"),
                  "ggplot")
  expect_error(multipatt_pq(data_fungi, fact = "Time"))
})
