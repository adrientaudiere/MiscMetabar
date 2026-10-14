skip_on_cran()

test_that("compare_pairs_pq() accepts fact and deprecates modality", {
  pq_lh <- subset_samples(data_fungi_mini, Height %in% c("Low", "High"))
  res_fact <- compare_pairs_pq(
    pq_lh,
    bifactor = "Height",
    fact = "Time",
    merge_sample_by = "Height"
  )
  lifecycle::expect_deprecated(
    res_mod <- compare_pairs_pq(
      pq_lh,
      bifactor = "Height",
      modality = "Time",
      merge_sample_by = "Height"
    )
  )
  expect_identical(res_fact, res_mod)
})

test_that("formattable_pq() accepts fact and deprecates modality", {
  skip_if_not_installed("formattable")
  res_fact <- formattable_pq(
    data_fungi_mini,
    fact = "Height",
    min_nb_seq_taxa = 5000,
    log10trans = TRUE
  )
  lifecycle::expect_deprecated(
    res_mod <- formattable_pq(
      data_fungi_mini,
      modality = "Height",
      min_nb_seq_taxa = 5000,
      log10trans = TRUE
    )
  )
  expect_s3_class(res_fact, "formattable")
  expect_identical(as.data.frame(res_fact), as.data.frame(res_mod))
})

test_that("hill_tuckey_pq() accepts fact and deprecates modality", {
  pq_sub <- prune_samples(
    sample_names(data_fungi_mini)[1:20],
    data_fungi_mini
  )
  res_fact <- hill_tuckey_pq(pq_sub, fact = "Height")
  lifecycle::expect_deprecated(
    res_mod <- hill_tuckey_pq(pq_sub, modality = "Height")
  )
  expect_s3_class(res_fact, "ggplot")
  expect_identical(res_fact$data, res_mod$data)
})

test_that("tax_datatable() accepts fact (column name or vector) and deprecates modality", {
  skip_on_os("windows")
  skip_if_not_installed("DT")
  pq_sub <- prune_taxa(taxa_names(data_fungi_mini)[1:10], data_fungi_mini)
  dt_col <- tax_datatable(pq_sub, fact = "Height")
  dt_vec <- tax_datatable(pq_sub, fact = pq_sub@sam_data$Height)
  lifecycle::expect_deprecated(
    dt_mod <- tax_datatable(pq_sub, modality = "Height")
  )
  expect_s3_class(dt_col, "datatables")
  expect_identical(dt_col$x$data, dt_vec$x$data)
  expect_identical(dt_col$x$data, dt_mod$x$data)
})

test_that("taxa_only_in_one_level() accepts fact and deprecates modality", {
  pq_no_na <- subset_samples(
    data_fungi_mini,
    !is.na(data_fungi_mini@sam_data$Height)
  )
  res_fact <- taxa_only_in_one_level(pq_no_na, fact = "Height", level = "High")
  lifecycle::expect_deprecated(
    res_mod <- taxa_only_in_one_level(
      pq_no_na,
      modality = "Height",
      level = "High"
    )
  )
  expect_identical(res_fact, res_mod)
})
