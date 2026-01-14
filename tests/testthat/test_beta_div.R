data(data_fungi)
data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High")) |>
  subset_samples(!is.na(Time) & !is.na(Height))

test_that("adonis_rarperm_pq works", {
  result <- adonis_rarperm_pq(data_subset, "Height", nperm = 9)
  expect_type(result, "list")
})

test_that("var_par_pq works", {
  result <- var_par_pq(data_subset,
    list_component = list(
      "Time" = c("Time"),
      "Size" = c("Height", "Diameter")
    ),
    dbrda_computation = TRUE
  )
  expect_s3_class(result, "varpart")
  expect_null(plot_var_part_pq(result))
})

test_that("var_par_rarperm_pq works", {
  result <- var_par_rarperm_pq(data_subset,
    list_component = list(
      "Time" = c("Time"),
      "Size" = c("Height", "Diameter")
    ), dbrda_computation = TRUE,
    nperm = 9
  )
  expect_type(result, "list")
  expect_s3_class(result, "varpart")
  expect_null(plot_var_part_pq(result,
    show_quantiles = TRUE,
    show_dbrda_signif = TRUE,
    show_dbrda_signif_pval = 0.1,
    id.size = 4,
    filter_quantile_zero = FALSE
  ))
})
