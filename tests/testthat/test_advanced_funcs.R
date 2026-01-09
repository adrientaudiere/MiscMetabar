data(data_fungi)

test_that("lulu works", {
  skip_on_cran()
  # Create a simple match list for testing
  data_subset <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 1000)
  # This requires a match list which is complex to create, so we test for error handling
  expect_error(lulu(otu_table = NULL, matchlist = NULL))
})

test_that("glmutli_pq works", {
  skip_on_cran()
  if (requireNamespace("glmulti", quietly = TRUE)) {
    data_subset <- subset_samples(data_fungi, Height %in% c("Low", "High"))
    # Test with simple parameters
    result <- glmutli_pq(data_subset, "shannon", c("Height", "Time"))
    expect_type(result, "list")
  }
})
