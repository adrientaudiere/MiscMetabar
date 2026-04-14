data("data_fungi")
data("GlobalPatterns", package = "phyloseq")
GP <- GlobalPatterns

data_fungi_2trees <-
  subset_samples(
    data_fungi,
    data_fungi@sam_data$Tree_name %in% c("A10-005", "AD30-abm-X")
  )
data_fungi_abun <-
  subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)

test_that("biplot_pq works", {
  skip_on_cran()
  expect_message(biplot_pq(data_fungi_2trees, merge_sample_by = "Tree_name"))
  expect_s3_class(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      plotly_version = TRUE
    ),
    "plotly"
  )
  expect_s3_class(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      geom_label = TRUE
    ),
    "ggplot"
  )
  expect_message(biplot_pq(
    data_fungi_2trees,
    fact = "Tree_name",
    merge_sample_by = "Tree_name",
    log10trans = FALSE,
    inverse_side = TRUE
  ))
  expect_error(
    biplot_pq(data_fungi, merge_sample_by = "Tree_name"),
    "biplot_pq needs only two samples"
  )
  expect_error(
    biplot_pq(data_fungi_2trees, fact = "Tree_name"),
    "biplot_pq needs only two samples"
  )
  expect_error(biplot_pq(data_fungi_2trees, merge_sample_by = "tRREE_name"))
  expect_s3_class(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      split_by_sample = TRUE
    ),
    "ggplot"
  )
  expect_s3_class(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      split_by_sample = TRUE,
      sample_border_col = "black",
      sample_border_width = 0.5
    ),
    "ggplot"
  )
  expect_error(biplot_pq(data_fungi_2trees, split_by_sample = TRUE))

  # color_rank
  expect_s3_class(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      color_rank = "Order"
    ),
    "ggplot"
  )
  expect_error(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      color_rank = "NotARank"
    ),
    "'color_rank' must be a column of tax_table"
  )

  # taxa_names_rank
  expect_s3_class(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      taxa_names_rank = "Genus"
    ),
    "ggplot"
  )
  expect_error(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      taxa_names_rank = "NotARank"
    ),
    "'taxa_names_rank' must be a column of tax_table"
  )

  # both together
  expect_s3_class(
    biplot_pq(
      data_fungi_2trees,
      merge_sample_by = "Tree_name",
      color_rank = "Order",
      taxa_names_rank = "Genus"
    ),
    "ggplot"
  )
})


test_that("multi_biplot_pq works with data_fungi dataset", {
  skip_on_cran()
  p1 <-
    multi_biplot_pq(data_fungi_abun, split_by = "Time", na_remove = FALSE)
  p2 <- multi_biplot_pq(data_fungi_abun, "Height")
  data_fungi_abun@sam_data$Random_pairs <-
    as.factor(sample(rep(1:(nsamples(data_fungi_abun) / 2), 2)))
  p3 <- multi_biplot_pq(data_fungi_abun, pairs = "Random_pairs")
  expect_s3_class(p1[[1]], "ggplot")
  expect_type(p1, "list")
  expect_s3_class(p2[[1]], "ggplot")
  expect_type(p2, "list")
  expect_s3_class(p3[[1]], "ggplot")
  expect_type(p3, "list")
  expect_length(p3, 85)
  expect_error(multi_biplot_pq(
    data_fungi_abun,
    pairs = "Random_pairs",
    split_by = "Time"
  ))
  expect_error(multi_biplot_pq(data_fungi_abun))
  expect_error(multi_biplot_pq(data_fungi_abun, pairs = "RandomPARR"))
  expect_error(multi_biplot_pq(data_fungi_abun, split_by = "TIMMEE"))
})
