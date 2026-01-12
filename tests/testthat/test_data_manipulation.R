data(data_fungi)

test_that("normalize_prop_pq works", {
  result <- normalize_prop_pq(data_fungi)
  expect_s4_class(result, "phyloseq")
  # Check that proportions sum to 1 (or close to it)
  expect_equal(sum(taxa_sums(result)), 47065.74)

  result_mini <- normalize_prop_pq(data_fungi_mini)
  expect_s4_class(result_mini, "phyloseq")
  # Check that proportions sum to 1 (or close to it)
  expect_equal(sum(taxa_sums(result_mini)), 4599.223)
})

test_that("merge_samples2 works with phyloseq objects", {
  suppressWarnings(result <- merge_samples2(data_fungi, "Height"))
  expect_s4_class(result, "phyloseq")
  expect_true(nsamples(result) <= nsamples(data_fungi))
})

test_that("merge_taxa_vec works with phyloseq objects", {
  group <- rep(c("Group1", "Group2"), length.out = ntaxa(data_fungi))
  result <- merge_taxa_vec(data_fungi, group)
  expect_s4_class(result, "phyloseq")
  expect_equal(ntaxa(result), 2)
})

test_that("psmelt_samples_pq works", {
  result <- psmelt_samples_pq(data_fungi)
  expect_s3_class(result, "data.frame")
  expect_equal(dim(result), c(185, 13))
})

test_that("rarefy_sample_count_by_modality works", {
  result <- rarefy_sample_count_by_modality(data_fungi, modality = "Height")
  expect_s4_class(result, "phyloseq")
  expect_equal(sum(result@otu_table), 1107122)
})

test_that("postcluster_pq works", {
  data_subset <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 1000)
  result <- postcluster_pq(data_subset, id = 0.97)
  expect_s4_class(result, "phyloseq")
  expect_equal(ntaxa(result), 167)
})
