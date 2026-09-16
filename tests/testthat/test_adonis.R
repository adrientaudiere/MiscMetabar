skip_on_cran()


data(enterotype, package = "phyloseq")

test_that("adonis function works fine", {
  expect_s3_class(adonis_pq(data_fungi, "Tree_name"), "anova")
  expect_s3_class(adonis_pq(data_fungi, "Height", na_remove = TRUE), "anova")
  expect_s3_class(
    adonis_pq(data_fungi, "Tree_name", correction_for_sample_size = TRUE),
    "anova"
  )
  expect_s3_class(
    adonis_pq(data_fungi, "Tree_name", rarefy_nb_seqs = TRUE),
    "anova"
  )
  expect_s3_class(
    adonis_pq(
      subset_samples(
        data_fungi,
        !is.na(data_fungi@sam_data$Time) & !is.na(data_fungi@sam_data$Height)
      ),
      "Time*Height"
    ),
    "anova"
  )
  expect_error(adonis_pq(data_fungi, "Time*Tree_name"))
  expect_error(adonis_pq(enterotype, "SeqTech*Enterotype"))
  expect_s3_class(
    adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE),
    "anova"
  )
  expect_error(adonis_pq(enterotype, "SecTech"))
  expect_error(adonis_pq(enterotype, "SeqTech", dist_method = "aitchison"))
})

test_that("adonis_pq forwards the by argument to vegan::adonis2()", {
  skip_on_cran()
  data(enterotype)
  res_terms <- adonis_pq(
    enterotype,
    "SeqTech+Enterotype",
    na_remove = TRUE,
    verbose = FALSE
  )
  res_model <- adonis_pq(
    enterotype,
    "SeqTech+Enterotype",
    na_remove = TRUE,
    by = NULL,
    verbose = FALSE
  )
  res_margin <- adonis_pq(
    enterotype,
    "SeqTech+Enterotype",
    na_remove = TRUE,
    by = "margin",
    verbose = FALSE
  )
  expect_true(all(c("SeqTech", "Enterotype") %in% rownames(res_terms)))
  expect_true(all(c("SeqTech", "Enterotype") %in% rownames(res_margin)))
  expect_true("Model" %in% rownames(res_model))
  expect_false(any(c("SeqTech", "Enterotype") %in% rownames(res_model)))
})
