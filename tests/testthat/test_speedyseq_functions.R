data(data_fungi)
data(data_fungi_mini)
data(enterotype)

test_that("merge_samples2 works with phyloseq objects", {
  # Create test data with merged samples
  ps <- enterotype
  sample_data(ps) <- sample_data(ps) %>%
    transform(Project.ClinicalStatus = Project:ClinicalStatus)
  expect_s4_class(merge_samples2(ps, "Project.ClinicalStatus"), "phyloseq")
  skip_on_cran()
  ps0 <- merge_samples2(ps, "Project.ClinicalStatus", fun_otu = mean)
  expect_s4_class(ps0, "phyloseq")
  expect_true(nsamples(ps0) < nsamples(ps))
})


test_that("merge_taxa_vec works with phyloseq objects", {
  if (requireNamespace("DECIPHER")) {
    # Create a simple grouping
    group_vec <- rep(1:3, length.out = ntaxa(data_fungi_mini))
    expect_s4_class(merge_taxa_vec(data_fungi_mini, group_vec), "phyloseq")
    skip_on_cran()
    result <- merge_taxa_vec(data_fungi_mini, group_vec)
    expect_true(ntaxa(result) <= 3)
    expect_true(ntaxa(result) < ntaxa(data_fungi_mini))
    
    # Test with NA in group
    group_vec_na <- group_vec
    group_vec_na[1] <- NA
    expect_warning(merge_taxa_vec(data_fungi_mini, group_vec_na))
  }
})


test_that("unique_or_na works correctly", {
  # Test with unique values
  expect_equal(unique_or_na(c("a", "a", "a")), "a")
  skip_on_cran()
  expect_equal(unique_or_na(c(1, 1, 1)), 1)
  
  # Test with non-unique values
  expect_true(is.na(unique_or_na(c("a", "b", "a"))))
  expect_true(is.na(unique_or_na(c(1, 2, 3))))
  
  # Test with factors
  f <- factor(c("a", "a", "b", "c"), ordered = TRUE)
  result <- unique_or_na(f)
  expect_true(is.na(result))
  
  result2 <- unique_or_na(f[1:2])
  expect_equal(as.character(result2), "a")
})


test_that("select_taxa works with otu_table", {
  otu <- otu_table(data_fungi_mini)
  taxa_to_select <- taxa_names(data_fungi_mini)[1:5]
  result <- select_taxa(otu, taxa_to_select)
  expect_s4_class(result, "otu_table")
  skip_on_cran()
  expect_equal(ntaxa(result), 5)
  
  # Test reorder
  result_no_reorder <- select_taxa(otu, taxa_to_select, reorder = FALSE)
  expect_s4_class(result_no_reorder, "otu_table")
})


test_that("select_taxa works with taxonomyTable", {
  tax <- tax_table(data_fungi_mini)
  taxa_to_select <- taxa_names(data_fungi_mini)[1:5]
  result <- select_taxa(tax, taxa_to_select)
  expect_s4_class(result, "taxonomyTable")
  skip_on_cran()
  expect_equal(ntaxa(result), 5)
})


test_that("select_taxa works with phyloseq objects", {
  taxa_to_select <- taxa_names(data_fungi_mini)[1:5]
  result <- select_taxa(data_fungi_mini, taxa_to_select)
  skip_on_cran()
  expect_true(is.null(result) || inherits(result, "phyloseq"))
})
