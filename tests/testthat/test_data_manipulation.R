data(data_fungi)
data(enterotype)

test_that("normalize_prop_pq works", {
  skip_on_cran()
  result <- normalize_prop_pq(data_fungi)
  expect_s4_class(result, "phyloseq")
  # Check that proportions sum to 1 (or close to it)
  sample_sums_result <- sample_sums(result)
  expect_true(all(abs(sample_sums_result - 1) < 0.01 | sample_sums_result == 0))
})

test_that("merge_samples2 works with phyloseq objects", {
  skip_on_cran()
  result <- merge_samples2(data_fungi, "Height")
  expect_s4_class(result, "phyloseq")
  expect_true(nsamples(result) <= nsamples(data_fungi))
})

test_that("merge_taxa_vec works with phyloseq objects", {
  skip_on_cran()
  # Create a simple grouping vector
  group <- rep(c("Group1", "Group2"), length.out = ntaxa(data_fungi))
  result <- merge_taxa_vec(data_fungi, group)
  expect_s4_class(result, "phyloseq")
  expect_true(ntaxa(result) <= ntaxa(data_fungi))
})

test_that("psmelt_samples_pq works", {
  skip_on_cran()
  result <- psmelt_samples_pq(data_fungi)
  expect_s3_class(result, "data.frame")
})

test_that("rarefy_sample_count_by_modality works", {
  skip_on_cran()
  result <- rarefy_sample_count_by_modality(data_fungi, modality = "Height")
  expect_s4_class(result, "phyloseq")
})

test_that("postcluster_pq works", {
  skip_on_cran()
  data_subset <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 1000)
  result <- postcluster_pq(data_subset, threshold = 0.97)
  expect_s4_class(result, "phyloseq")
})
