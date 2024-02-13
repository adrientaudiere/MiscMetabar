data("GlobalPatterns", package = "phyloseq")
data("enterotype", package = "phyloseq")

GP <- GlobalPatterns
data_fungi_2trees <-
  subset_samples(
    data_fungi,
    data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X")
  )
GP_archae <-
  subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")

test_that("rotl_pq works with data_fungi dataset", {
  skip_on_os("windows")
  skip_on_cran()
  library("rotl")
  expect_s3_class(suppressWarnings(tr <-
    rotl_pq(data_fungi, species_colnames = "Genus_species")), "phylo")
  expect_s3_class(suppressWarnings(
    rotl_pq(
      data_fungi,
      species_colnames = "Genus_species",
      context_name = "Ascomycetes"
    )
  ), "phylo")
  expect_silent(plot(tr))
})

test_that("heat_tree_pq works with data_fungi dataset", {
  skip_on_cran()
  library(metacoder)
  expect_silent(suppressMessages(ht <- heat_tree_pq(data_fungi_mini)))
  expect_s3_class(ht, "ggplot")
  expect_s3_class(
    heat_tree_pq(data_fungi_mini, taxonomic_level = 1:4),
    "ggplot"
  )
})

GPsubset <- subset_taxa(
  GlobalPatterns,
  GlobalPatterns@tax_table[, 1] == "Bacteria"
)
test_that("heat_tree_pq works with GlobalPatterns dataset", {
  skip_on_cran()
  library(metacoder)
  expect_silent(suppressMessages(ht <- heat_tree_pq(GPsubset)))
  expect_silent(suppressMessages(
    ht <-
      heat_tree_pq(
        GPsubset,
        node_size = n_obs,
        node_color = n_obs,
        node_label = taxon_names,
        tree_label = taxon_names,
        node_size_trans = "log10 area"
      )
  ))
  expect_s3_class(ht, "ggplot")
})


test_that("plot_tax_pq works with data_fungi dataset", {
  expect_silent(suppressMessages(
    pt <-
      plot_tax_pq(
        data_fungi_mini,
        "Time",
        merge_sample_by = "Time",
        taxa_fill = "Class",
        add_info = FALSE
      )
  ))
  skip_on_cran()
  expect_silent(suppressMessages(
    pt <-
      plot_tax_pq(
        data_fungi_mini,
        "Time",
        merge_sample_by = "Time",
        taxa_fill = "Class",
        nb_print_value = 2
      )
  ))
  expect_s3_class(pt, "ggplot")
  expect_silent(suppressMessages(
    pt <-
      plot_tax_pq(data_fungi_mini, "Time", taxa_fill = "Class")
  ))
  expect_s3_class(pt, "ggplot")
  expect_silent(suppressMessages(
    pt <-
      plot_tax_pq(
        data_fungi_mini,
        "Time",
        merge_sample_by = "Time",
        taxa_fill = "Class",
        type = "nb_asv",
        add_info = FALSE
      )
  ))
  expect_silent(suppressMessages(
    pt <-
      plot_tax_pq(
        data_fungi_mini,
        "Time",
        merge_sample_by = "Time",
        taxa_fill = "Class",
        type = "nb_asv"
      )
  ))
  expect_s3_class(pt, "ggplot")
  expect_silent(suppressMessages(
    pt <-
      plot_tax_pq(
        data_fungi_mini,
        "Time",
        merge_sample_by = "Time",
        taxa_fill = "Class",
        type = "both",
        add_info = FALSE
      )
  ))
  expect_silent(suppressMessages(
    pt <-
      plot_tax_pq(
        data_fungi_mini,
        "Time",
        merge_sample_by = "Time",
        taxa_fill = "Class",
        type = "both"
      )
  ))
  expect_s3_class(pt[[1]], "ggplot")
  expect_silent(suppressMessages(
    plot_tax_pq(
      data_fungi_mini,
      "Time",
      merge_sample_by = "Time",
      taxa_fill = "Class",
      na_remove = TRUE
    )
  ))
  expect_silent(suppressMessages(
    plot_tax_pq(
      data_fungi_mini,
      "Time",
      merge_sample_by = "Time",
      taxa_fill = "Order",
      clean_pq = FALSE
    )
  ))
  expect_silent(suppressMessages(
    plot_tax_pq(
      data_fungi_mini,
      "Height",
      merge_sample_by = "Height",
      taxa_fill = "Order"
    )
  ))
})


test_that("multitax_bar_pq works with data_fungi_sp_known dataset", {
  expect_s3_class(
    multitax_bar_pq(data_fungi_mini, "Phylum", "Class", "Order", "Time"),
    "ggplot"
  )
  skip_on_cran()
  expect_s3_class(
    multitax_bar_pq(data_fungi_mini, "Phylum", "Class", "Order"),
    "ggplot"
  )
  expect_s3_class(
    multitax_bar_pq(
      data_fungi_mini,
      "Phylum",
      "Class",
      "Order",
      nb_seq = FALSE,
      log10trans = FALSE
    ),
    "ggplot"
  )
  expect_s3_class(
    multitax_bar_pq(data_fungi_mini,
      "Phylum",
      "Class",
      "Order",
      log10trans = FALSE
    ),
    "ggplot"
  )
  expect_error(print(
    multitax_bar_pq(data_fungi_mini, "Class", "Genus", "Order")
  ))
  expect_error(print(
    multitax_bar_pq(data_fungi_mini, "Phylum", "Class", "Order", "TIMESS")
  ))
})


test_that("multitax_bar_pq works with GlobalPatterns dataset", {
  expect_s3_class(
    multitax_bar_pq(GP_archae, "Phylum", "Class", "Order", "SampleType"),
    "ggplot"
  )
  skip_on_cran()
  expect_s3_class(
    multitax_bar_pq(GP_archae, "Phylum", "Class", "Order", nb_seq = FALSE),
    "ggplot"
  )
  expect_s3_class(
    multitax_bar_pq(
      GP_archae,
      "Phylum",
      "Class",
      "Order",
      nb_seq = FALSE,
      log10trans = FALSE
    ),
    "ggplot"
  )
  expect_s3_class(
    multitax_bar_pq(GP_archae,
      "Phylum",
      "Class",
      "Order",
      log10trans = FALSE
    ),
    "ggplot"
  )
  expect_error(print(multitax_bar_pq(GP_archae, "Class", "Genus", "Order")))
  expect_error(print(
    multitax_bar_pq(GP_archae, "Phylum", "Class", "Order", "UNKOWNS")
  ))
})


test_that("rigdes_pq work with data_fungi dataset", {
  expect_s3_class(
    ridges_pq(data_fungi_mini,
      "Time",
      alpha = 0.5,
      log10trans = FALSE
    ) + xlim(c(0, 1000)),
    "ggplot"
  )
  skip_on_cran()
  expect_s3_class(
    ridges_pq(data_fungi_mini, "Time", alpha = 0.5),
    "ggplot"
  )
  expect_s3_class(
    ridges_pq(
      data_fungi_mini,
      "Time",
      nb_seq = FALSE,
      log10trans = FALSE
    ),
    "ggplot"
  )
  expect_s3_class(
    ridges_pq(
      clean_pq(
        subset_taxa(data_fungi_sp_known, Phylum == "Basidiomycota")
      ),
      "Time"
    ),
    "ggplot"
  )
  expect_s3_class(
    ridges_pq(clean_pq(
      subset_taxa(data_fungi_sp_known, Phylum == "Basidiomycota")
    ), "Time", alpha = 0.6, scale = 0.9),
    "ggplot"
  )
  expect_s3_class(
    ridges_pq(
      clean_pq(subset_taxa(
        data_fungi_sp_known, Phylum == "Basidiomycota"
      )),
      "Time",
      jittered_points = TRUE,
      position = ggridges::position_points_jitter(width = 0.05, height = 0),
      point_shape = "|",
      point_size = 3,
      point_alpha = 1,
      alpha = 0.7,
      scale = 0.8
    ),
    "ggplot"
  )
  expect_error(
    ridges_pq(clean_pq(
      subset_taxa(data_fungi_sp_known, Phylum == "Basidiomycota")
    ))
  )
})



test_that("treemap_pq work with data_fungi_sp_known dataset", {
  expect_s3_class(
    treemap_pq(
      clean_pq(
        data_fungi_mini
      ),
      "Order", "Class",
      plot_legend = TRUE
    ),
    "ggplot"
  )
  skip_on_cran()
  expect_s3_class(
    treemap_pq(
      clean_pq(
        data_fungi_mini
      ),
      "Order", "Class",
      log10trans = FALSE
    ),
    "ggplot"
  )
  expect_s3_class(
    treemap_pq(
      data_fungi_mini,
      "Order",
      "Class",
      nb_seq = FALSE,
      log10trans = FALSE
    ),
    "ggplot"
  )
  expect_s3_class(
    treemap_pq(
      data_fungi_mini,
      "Order",
      "Class",
      nb_seq = FALSE,
      log10trans = TRUE
    ),
    "ggplot"
  )
})

test_that("tax_bar_pq work with data_fungi dataset", {
  expect_s3_class(tax_bar_pq(data_fungi_mini, taxa = "Class"), "ggplot")
  skip_on_cran()
  expect_s3_class(tax_bar_pq(data_fungi_mini, taxa = "Class", fact = "Time"), "ggplot")
  expect_s3_class(
    tax_bar_pq(
      data_fungi_mini,
      taxa = "Class",
      fact = "Time",
      nb_seq = FALSE
    ),
    "ggplot"
  )
  expect_s3_class(
    tax_bar_pq(
      data_fungi_mini,
      taxa = "Class",
      fact = "Time",
      nb_seq = FALSE,
      percent_bar = TRUE
    ),
    "ggplot"
  )
})

test_that("add_funguild_info and plot_guild_pq work with data_fungi_mini dataset", {
  skip_on_cran()
  expect_s4_class(
    df <-
      add_funguild_info(
        subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 5000),
        taxLevels = c(
          "Domain",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species"
        )
      ),
    "phyloseq"
  )
  expect_error(df <-
    add_funguild_info(
      subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 5000),
      taxLevels = c(
        "PHYLLUUM",
        "Phylum",
        "Class",
        "Order",
        "Family"
      )
    ))
  expect_s3_class(plot_guild_pq(df, clean_pq = TRUE), "ggplot")
  expect_s3_class(plot_guild_pq(df, clean_pq = FALSE), "ggplot")
})


test_that("build_phytree_pq work with data_fungi dataset", {
  skip_on_os("windows")
  skip_on_cran()
  df <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 19000)
  expect_type(df_tree <- build_phytree_pq(df, nb_bootstrap = 2, rearrangement = "stochastic"), "list")
  expect_type(df_tree <- build_phytree_pq(df, nb_bootstrap = 2, rearrangement = "ratchet"), "list")
  expect_error(build_phytree_pq(df, nb_bootstrap = 2, rearrangement = "PRAtchet"))
  expect_error(build_phytree_pq(GP, nb_bootstrap = 2))

  expect_type(df_tree <- build_phytree_pq(df, nb_bootstrap = 2), "list")
  expect_equal(length(df_tree), 6)
  expect_type(df_tree_wo_bootstrap <- build_phytree_pq(df, nb_bootstrap = 0), "list")
  expect_equal(length(df_tree_wo_bootstrap), 3)
  expect_s3_class(df_tree$NJ, "phylo")
  expect_s3_class(df_tree$UPGMA, "phylo")
  expect_s3_class(df_tree$ML, "pml")
  expect_s3_class(df_tree$NJ_bs, "multiPhylo")
  expect_s3_class(df_tree$UPGMA_bs, "multiPhylo")
  expect_s3_class(df_tree$ML_bs, "multiPhylo")
})
