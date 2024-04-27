skip_on_cran()

data(data_fungi)
data(enterotype, package = "phyloseq")

test_that("adonis function works fine", {
  expect_s3_class(adonis_pq(data_fungi, "Tree_name"), "anova")
  expect_s3_class(adonis_pq(data_fungi, "Height", na_remove = TRUE), "anova")
  expect_s3_class(adonis_pq(data_fungi, "Tree_name", correction_for_sample_size = TRUE), "anova")
  expect_s3_class(adonis_pq(data_fungi, "Tree_name", rarefy_nb_seqs = TRUE), "anova")
  expect_s3_class(adonis_pq(subset_samples(data_fungi, !is.na(data_fungi@sam_data$Time)), "Time*Tree_name"), "anova")
  expect_error(adonis_pq(data_fungi, "Time*Tree_name"))
  expect_error(adonis_pq(enterotype, "SeqTech*Enterotype"))
  expect_s3_class(adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE), "anova")
  expect_error(adonis_pq(enterotype, "SecTech"))
  expect_error(adonis_pq(enterotype, "SeqTech", dist_method = "aitchison"))
})
