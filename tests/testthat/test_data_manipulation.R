data(data_fungi)

test_that("normalize_prop_pq works", {
  result <- normalize_prop_pq(data_fungi)
  expect_s4_class(result, "phyloseq")
  # Check that proportions sum to 1 (or close to it)
  expect_equal(sum(taxa_sums(result)), 47065, tolerance = 1)

  result_mini <- normalize_prop_pq(data_fungi_mini)
  expect_s4_class(result_mini, "phyloseq")
  # Check that proportions sum to 1 (or close to it)
  expect_equal(sum(taxa_sums(result_mini)), 4599, tolerance = 1)
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
  expect_identical(ntaxa(result), 2L)
})

test_that("merge_taxa_vec handles NA groups with a warning", {
  group <- rep(c("G1", "G2"), length.out = ntaxa(data_fungi))
  group[1] <- NA
  expect_warning(merge_taxa_vec(data_fungi, group), "missing values")
})

test_that("merge_taxa_vec with tax_adjust=2 works", {
  group <- rep(c("G1", "G2"), length.out = ntaxa(data_fungi))
  result <- merge_taxa_vec(data_fungi, group, tax_adjust = 2L)
  expect_s4_class(result, "phyloseq")
  expect_identical(ntaxa(result), 2L)
})

test_that("merge_taxa_vec with rank_propagation=FALSE works", {
  group <- rep(c("G1", "G2"), length.out = ntaxa(data_fungi))
  result <- merge_taxa_vec(data_fungi, group, rank_propagation = FALSE)
  expect_s4_class(result, "phyloseq")
  expect_identical(ntaxa(result), 2L)
})

test_that("merge_samples2 handles NA group values with a warning", {
  grp <- data_fungi@sam_data[["Height"]]
  grp[1] <- NA
  expect_warning(merge_samples2(data_fungi, grp), "missing values")
})

test_that("merge_samples2 works with non-sum fun_otu", {
  result <- suppressWarnings(
    merge_samples2(data_fungi, "Height", fun_otu = mean)
  )
  expect_s4_class(result, "phyloseq")
  expect_true(nsamples(result) <= nsamples(data_fungi))
})

test_that("psmelt_samples_pq works", {
  result <- psmelt_samples_pq(data_fungi)
  expect_s3_class(result, "data.frame")
  expect_identical(dim(result), c(185L, 13L))
})

test_that("rarefy_sample_count_by_modality works", {
  result <- rarefy_sample_count_by_modality(
    data_fungi,
    fact = "Height",
    rngseed = 42
  )
  expect_s4_class(result, "phyloseq")
  expect_equal(sum(result@otu_table), 1065385)
})

data_subset <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 1000)

test_that("postcluster_pq works with method ='clusterize'", {
  set.seed(42)
  result <- postcluster_pq(data_subset, id = 0.97)
  expect_s4_class(result, "phyloseq")
  expect_identical(ntaxa(result), 167L)
})


test_that("postcluster_pq works with method ='vsearch'", {
  skip_on_cran()
  set.seed(42)
  result <- postcluster_pq(data_subset, method = "vsearch")
  expect_s4_class(result, "phyloseq")
  expect_identical(ntaxa(result), 168L)
})


test_that("postcluster_pq works with method ='vsearch' and
   vsearch_cluster_method='cluster_fast'", {
  skip_on_cran()
  set.seed(42)
  result <- postcluster_pq(
    data_subset,
    method = "vsearch",
    vsearch_cluster_method = "--cluster_fast"
  )
  expect_s4_class(result, "phyloseq")
  expect_identical(ntaxa(result), 170L)
})

test_that("postcluster_pq works with method ='swarm'", {
  skip_on_cran()
  set.seed(42)
  result <- postcluster_pq(data_subset, method = "swarm")
  expect_s4_class(result, "phyloseq")
  expect_identical(ntaxa(result), 216L)
})
