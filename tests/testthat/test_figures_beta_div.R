data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns
data_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
data_basidio_2trees <-
  subset_samples(data_basidio, Tree_name %in% c("A10-005", "AD30-abm-X"))
GP_archae <-
  subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")

data_fungi_woNA4time <- subset_samples(data_fungi_mini, !is.na(Time))
res_mt <-
  mt(
    data_fungi_woNA4time,
    "Time",
    method = "fdr",
    test = "f",
    B = 300
  )

test_that("circle_pq works", {
  skip_on_cran()
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
  if (requireNamespace("phyloseqGraphTest")) {
    expect_silent(graph_test_pq(data_fungi_mini, fact = "Tree_name"))
    skip_on_cran()
    expect_silent(graph_test_pq(data_fungi_mini, fact = "Tree_name", na_remove = TRUE))
    expect_silent(graph_test_pq(data_fungi_mini, fact = "Tree_name", return_plot = FALSE))
    expect_silent(graph_test_pq(
      subset_samples(data_fungi_mini, !is.na(data_fungi_mini@sam_data$Time)),
      fact = "Time",
      merge_sample_by = "Tree_name"
    ))
    expect_error(graph_test_pq(data_fungi_mini, fact = "Height"))
    expect_message(graph_test_pq(data_fungi_mini, fact = "Height", na_remove = TRUE))
    expect_error(graph_test_pq(enterotype, fact = "Enterotype"))
    expect_error(graph_test_pq(data_fungi_mini, fact = "tRREE_name"))
  }
})


test_that("plot_mt works", {
  expect_s3_class(res_mt, "data.frame")
  expect_s3_class(suppressWarnings(plot_mt(res_mt)), "ggplot")
  skip_on_cran()
  expect_s3_class(suppressWarnings(plot_mt(
    res_mt,
    taxa = "Genus", color_tax = "Order"
  )), "ggplot")
})

test_that("sankey_pq works with GlobalPatterns dataset", {
  if (requireNamespace("networkD3")) {
    expect_silent(sankey_pq(GP_archae))
    skip_on_cran()
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
  }
})

test_that("sankey_pq works with data_fungi_mini dataset", {
  if (requireNamespace("networkD3")) {
    expect_silent(sankey_pq(data_fungi_mini))
    skip_on_cran()
    expect_s3_class(sankey_pq(data_fungi_mini), "htmlwidget")
    expect_s3_class(sankey_pq(data_fungi_mini), "sankeyNetwork")
    expect_silent(suppressWarnings(sankey_pq(data_fungi_mini, fact = "Height")))
    expect_silent(sankey_pq(
      data_fungi_mini,
      taxa = 3:7,
      min_prop_tax = 0.01,
      units = "sequences"
    ))
    expect_silent(sankey_pq(
      data_fungi_mini,
      taxa = 1:4,
      min_prop_tax = 0.01,
      add_nb_seq = TRUE
    ))
    expect_silent(
      sankey_pq(
        data_fungi_mini,
        taxa = 1:4,
        min_prop_tax = 0.001,
        add_nb_seq = TRUE,
        tax2remove = "Undefined"
      )
    )
    expect_silent(
      sankey_pq(
        data_fungi_mini,
        taxa = 1:4,
        add_nb_seq = TRUE,
        units = "sequences",
        symbol2sub = NULL
      )
    )
    expect_warning(
      sankey_pq(
        data_fungi_mini,
        taxa = 1:4,
        min_prop_tax = 0.01,
        add_nb_seq = TRUE,
        units = "sequences",
        symbol2sub = NA
      )
    )
    expect_error(sankey_pq(data_fungi_mini, "HEIGHT"))
  }
})

test_that("venn_pq works with data_fungi_mini dataset", {
  skip_on_os("windows")
  skip_on_cran()
  library("grid")
  expect_silent(venn_pq(data_fungi_mini, "Height"))
  expect_silent(suppressMessages(venn_pq(
    clean_pq(
      data_fungi_mini,
      force_taxa_as_rows = TRUE,
      silent = TRUE
    ),
    "Height"
  )))
  expect_silent(venn_pq(data_fungi_mini, "Height", min_nb_seq = 10))
  expect_silent(venn_pq(data_fungi_mini, "Height", print_values = FALSE))
  expect_silent(venn_pq(data_fungi_mini, "Height", print_values = FALSE) +
    scale_fill_hue())
  expect_silent(venn_pq(data_fungi_mini, "Height", print_values = TRUE) +
    scale_fill_hue())
  expect_error(venn_pq(data_fungi_mini))
  expect_error(venn_pq(data_fungi_mini@otu_table, "Height"))
  expect_type(
    venn_pq(data_fungi_mini, "Height", print_values = TRUE) + scale_fill_hue(),
    "NULL"
  )
  expect_s3_class(
    venn_pq(data_fungi_mini, "Height", print_values = FALSE) + scale_fill_hue(),
    "ggplot"
  )
})

test_that("ggvenn_pq works with data_fungi_mini dataset", {
  if (requireNamespace("ggVennDiagram")) {
    expect_message(ggvenn_pq(data_fungi_mini, "Height"))
    skip_on_cran()
    expect_message(ggvenn_pq(data_fungi_mini, "Height", rarefy_before_merging = TRUE))
    expect_message(suppressWarnings(ggvenn_pq(data_fungi_mini, "Height", rarefy_after_merging = TRUE)))
    expect_message(ggvenn_pq(data_fungi_mini, "Height", add_nb_seq = TRUE))
    expect_silent(suppressMessages(ggvenn_pq(data_fungi_mini, "Height", rarefy_nb_seqs = TRUE)))
    expect_message(ggvenn_pq(data_fungi_mini, "Height", min_nb_seq = 2))
    expect_message(ggvenn_pq(data_fungi_mini, "Height", taxonomic_rank = 4))
    expect_silent(suppressMessages(ggvenn_pq(data_fungi_mini, "Height", split_by = "Time")))
    expect_error(ggvenn_pq(data_fungi_mini))
    expect_s3_class(suppressMessages(ggvenn_pq(data_fungi_mini, "Height")), "ggplot")
    expect_error(ggvenn_pq(data_fungi_mini@otu_table, "Height"))
  }
})


test_that("upset_pq works with data_fungi dataset", {
if (requireNamespace("tidyr") && requireNamespace("ComplexUpset") && packageVersion("ggplot2") < "4.0.0") {
    expect_silent(suppressMessages(upset_pq(data_fungi_mini, "Height")))
    skip_on_cran()
    expect_s3_class(upset_pq(data_fungi_mini, "Height", taxa_fill = "Class"), "ggplot")
    expect_s3_class(upset_pq(data_fungi_mini, "Height"), "ggplot")
    expect_s3_class(
      upset_pq(
        data_fungi_mini,
        "Height",
        na_remove = TRUE,
        rarefy_after_merging = TRUE
      ),
      "ggplot"
    )

    expect_s3_class(upset_pq(data_fungi_mini, "Time"), "ggplot")
    expect_s3_class(upset_pq(data_fungi_mini, "Time", min_nb_seq = 10), "ggplot")
    expect_s3_class(
      upset_pq(data_fungi_mini, "Time",
        numeric_fonction = mean,
        na_remove = FALSE
      ),
      "ggplot"
    )

    expect_error(upset_pq(data_fungi_mini))
  }
})

test_that("upset_test_pq works with data_fungi_mini dataset", {
  if (requireNamespace("tidyr") && requireNamespace("ComplexUpset")) {
    expect_s3_class(upset_test_pq(data_fungi_mini, "Height"), "data.frame")
    skip_on_cran()
    expect_s3_class(upset_test_pq(data_fungi_mini, "Time"), "data.frame")
    expect_s3_class(
      upset_test_pq(data_fungi_mini, "Time", min_nb_seq = 10),
      "data.frame"
    )
    expect_s3_class(
      upset_test_pq(data_fungi_mini, "Time", numeric_fonction = mean),
      "data.frame"
    )
    expect_s3_class(
      upset_test_pq(
        data_fungi_mini,
        "Time",
        numeric_fonction = mean,
        var_to_test = c("OTU", "Guild", "Genus")
      ),
      "data.frame"
    )
    expect_error(upset_test_pq(data_fungi_mini, "Height", var_to_test = "GUILDDDS"))
    expect_error(upset_test_pq(data_fungi_mini))
  }
})

test_that("plot_LCBD_pq works with data_fungi dataset", {
  skip_on_cran()
  expect_s3_class(
    suppressWarnings(plot_LCBD_pq(
      data_fungi_mini,
      nperm = 100,
      only_plot_significant = FALSE
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(plot_LCBD_pq(
      data_fungi_mini,
      nperm = 100,
      only_plot_significant = TRUE,
      pval = 0.2
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(plot_LCBD_pq(
      data_fungi_mini,
      nperm = 100,
      only_plot_significant = TRUE,
      p_adjust_method = "holm",
      sam_variables = c("Time", "Height")
    )),
    "ggplot"
  )
})


test_that("LCBD_pq works with data_fungi_mini dataset", {
  skip_on_cran()
  expect_s3_class(LCBD_pq(data_fungi_mini, nperm = 100), "beta.div")
  expect_s3_class(
    LCBD_pq(data_fungi_mini, nperm = 100, method = "jaccard"),
    "beta.div"
  )
})

test_that("plot_LCBD_pq works with data_fungi_mini dataset", {
  skip_on_cran()
  expect_s3_class(
    suppressWarnings(plot_LCBD_pq(
      data_fungi_mini,
      nperm = 100,
      only_plot_significant = FALSE
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(plot_LCBD_pq(
      data_fungi_mini,
      nperm = 100,
      only_plot_significant = TRUE,
      pval = 0.2
    )),
    "ggplot"
  )
  expect_s3_class(
    suppressWarnings(plot_LCBD_pq(
      data_fungi_mini,
      nperm = 100,
      only_plot_significant = TRUE,
      p_adjust_method = "holm",
      sam_variables = c("Time", "Height")
    )),
    "ggplot"
  )
})


test_that("plot_SCBD_pq works with data_fungi_mini dataset", {
  skip_on_cran()
  expect_s3_class(suppressWarnings(plot_SCBD_pq(data_fungi_mini)), "ggplot")
  expect_s3_class(
    suppressWarnings(plot_SCBD_pq(
      data_fungi_mini,
      tax_level = "Class",
      tax_col = "Phylum",
      min_SCBD = 0
    )),
    "ggplot"
  )
})

test_that("multipatt_pq works with data_fungi_mini dataset", {
  skip_on_os("windows")
  skip_on_cran()
  expect_s3_class(
    multipatt_pq(subset_samples(data_fungi_mini, !is.na(Time)),
      fact = "Time"
    ),
    "ggplot"
  )
  expect_error(multipatt_pq(data_fungi_mini, fact = "Time"))
})

test_that("multipatt_pq works with data_fungi_mini dataset", {
  skip_on_os("windows")
  skip_on_cran()
  expect_type(suppressMessages(suppressWarnings(res_height <- ancombc_pq(
    subset_taxa_pq(
      data_fungi_sp_known,
      taxa_sums(data_fungi_sp_known) > 5000
    ),
    fact = "Height",
    levels_fact = c("Low", "High"),
    verbose = TRUE
  ))), "list")
  expect_s3_class(res_height$bias_correct_log_table, "data.frame")
  expect_equal(dim(res_height$res), c(6, 17))

  expect_type(suppressMessages(suppressWarnings(res_time <- ancombc_pq(
    subset_taxa_pq(
      data_fungi_sp_known,
      taxa_sums(data_fungi_sp_known) > 5000
    ),
    fact = "Time",
    levels_fact = c("0", "15"),
    tax_level = "Family",
    verbose = TRUE
  ))), "list")

  expect_s3_class(res_time$ss_tab, "data.frame")
  expect_equal(dim(res_time$res), c(12, 17))
})
