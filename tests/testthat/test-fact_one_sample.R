data("GlobalPatterns", package = "phyloseq")

# Create a single-level factor phyloseq: subset to 1 sample type
pq_one <- subset_samples(
  GlobalPatterns,
  SampleType == "Soil"
)
pq_one@sam_data[["SingleFact"]] <- factor("GroupA")

# Category A: functions that SHOULD WORK with 1 factor level

test_that("tax_bar_pq works with 1-level factor", {
  skip_on_cran()
  expect_s3_class(
    tax_bar_pq(pq_one, fact = "SingleFact", taxa = "Phylum"),
    "ggplot"
  )
})

test_that("plot_tax_pq works with 1-level factor", {
  skip_on_cran()
  expect_s3_class(
    suppressMessages(
      plot_tax_pq(pq_one, fact = "SingleFact", taxa_fill = "Phylum")
    ),
    "ggplot"
  )
})

test_that("ridges_pq works with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("ggridges", quietly = TRUE)) {
    expect_s3_class(
      suppressWarnings(
        ridges_pq(pq_one, fact = "SingleFact", tax_level = "Phylum")
      ),
      "ggplot"
    )
  }
})

test_that("multitax_bar_pq works with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("ggh4x", quietly = TRUE)) {
    expect_s3_class(
      multitax_bar_pq(
        pq_one,
        lvl1 = "Phylum",
        lvl2 = "Class",
        lvl3 = "Order",
        fact = "SingleFact"
      ),
      "ggplot"
    )
  }
})

test_that("sankey_pq works with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("networkD3", quietly = TRUE)) {
    expect_no_error(
      suppressMessages(sankey_pq(pq_one, fact = "SingleFact", taxa = 1:3))
    )
  }
})

test_that("circle_pq works with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("circlize", quietly = TRUE)) {
    expect_no_error(
      suppressMessages(circle_pq(pq_one, fact = "SingleFact", nproc = 1))
    )
  }
})

test_that("ggaluv_pq works with single-sample-type phyloseq", {
  skip_on_cran()
  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    expect_s3_class(
      suppressMessages(ggaluv_pq(pq_one)),
      "ggplot"
    )
  }
})

test_that("distri_1_taxa works with 1-level factor", {
  skip_on_cran()
  taxa_name <- phyloseq::taxa_names(pq_one)[1]
  expect_no_error(
    res <- distri_1_taxa(pq_one, "SingleFact", taxa_name)
  )
  expect_s3_class(res, "data.frame")
})

test_that("are_modality_even_depth works with 1-level factor", {
  skip_on_cran()
  expect_no_error(
    res <- are_modality_even_depth(pq_one, "SingleFact")
  )
})

test_that("rarefy_sample_count_by_modality works with 1-level factor", {
  skip_on_cran()
  expect_no_error(
    suppressMessages(
      res <- rarefy_sample_count_by_modality(
        pq_one,
        "SingleFact",
        rngseed = 1
      )
    )
  )
  expect_s4_class(res, "phyloseq")
})

# Category B: functions that SHOULD GIVE INFORMATIVE ERROR with 1 level

test_that("hill_pq errors informatively with 1-level factor", {
  skip_on_cran()
  expect_error(
    hill_pq(pq_one, fact = "SingleFact"),
    "at least (two|2)"
  )
})

test_that("hill_test_rarperm_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("ggstatsplot", quietly = TRUE)) {
    expect_error(
      hill_test_rarperm_pq(pq_one, fact = "SingleFact", nperm = 1),
      "at least (two|2)"
    )
  }
})

test_that("graph_test_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("phyloseqGraphTest", quietly = TRUE)) {
    expect_error(
      graph_test_pq(pq_one, fact = "SingleFact"),
      "at least (two|2)"
    )
  }
})

test_that("multipatt_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("indicspecies", quietly = TRUE)) {
    expect_error(
      multipatt_pq(pq_one, fact = "SingleFact"),
      "at least (two|2)"
    )
  }
})

test_that("ancombc_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("ANCOMBC", quietly = TRUE)) {
    expect_error(
      ancombc_pq(pq_one, fact = "SingleFact"),
      "at least (two|2)"
    )
  }
})

test_that("ggbetween_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("ggstatsplot", quietly = TRUE)) {
    expect_error(
      ggbetween_pq(pq_one, fact = "SingleFact"),
      "at least (two|2)"
    )
  }
})

test_that("venn_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("venneuler", quietly = TRUE)) {
    expect_error(
      venn_pq(pq_one, fact = "SingleFact"),
      "at least (two|2)"
    )
  }
})

test_that("ggvenn_pq errors informatively with 1-level factor", {
  skip_on_cran()
  expect_error(
    ggvenn_pq(pq_one, fact = "SingleFact"),
    "at least (two|2)"
  )
})

test_that("upset_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    expect_error(
      upset_pq(pq_one, fact = "SingleFact"),
      "at least (two|2)"
    )
  }
})

test_that("biplot_pq errors informatively with 1-level factor", {
  skip_on_cran()
  expect_error(
    biplot_pq(pq_one, fact = "SingleFact", merge_sample_by = "SingleFact"),
    "(two|2)"
  )
})

test_that("accu_plot errors informatively with 1-level factor", {
  skip_on_cran()
  expect_error(
    accu_plot(pq_one, fact = "SingleFact"),
    "at least (two|2)"
  )
})

test_that("accu_plot_balanced_modality errors with 1-level factor", {
  skip_on_cran()
  expect_error(
    accu_plot_balanced_modality(pq_one, fact = "SingleFact"),
    "at least (two|2)"
  )
})

test_that("plot_tsne_pq errors informatively with 1-level factor", {
  skip_on_cran()
  if (requireNamespace("Rtsne", quietly = TRUE)) {
    expect_error(
      plot_tsne_pq(pq_one, fact = "SingleFact"),
      "at least (two|2)"
    )
  }
})
