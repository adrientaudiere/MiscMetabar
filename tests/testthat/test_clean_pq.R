data("data_fungi")

data_fungi_test <- data_fungi
data_fungi_test@otu_table[, 1] <- rep(0, nrow(data_fungi_test@otu_table))
data_fungi_test@otu_table[10, ] <- rep(0, ncol(data_fungi_test@otu_table))

test_that("clean_pq return a phyloseq object after cleaning empty taxa and samples", {
  expect_s4_class(data_fungi_test, "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test), "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test, verbose = TRUE), "phyloseq")
})

test_that("clean_pq clean empty taxa and samples", {
  expect_message(clean_pq(data_fungi_test), "Cleaning suppress 2 taxa and 1 samples")
  expect_equal(nrow(data_fungi_test@otu_table) -
    nrow(clean_pq(data_fungi_test)@otu_table), 1)
  expect_equal(ncol(data_fungi_test@otu_table) -
    ncol(clean_pq(data_fungi_test)@otu_table), 2)
})

test_that("clean_pq force taxa in column", {
  expect_true(taxa_are_rows(clean_pq(data_fungi, force_taxa_as_rows = TRUE)))
  expect_false(taxa_are_rows(clean_pq(data_fungi, force_taxa_as_columns = TRUE)))
})

data_fungi_test2 <- data_fungi_test
taxa_names(data_fungi_test2) <- paste0("ASV", 1:ntaxa(data_fungi_test2))
test_that("clean_pq works fine with bad taxa_names", {
  expect_s4_class(clean_pq(data_fungi_test2), "phyloseq")
    expect_message(clean_pq(data_fungi_test2), "Cleaning suppress 2 taxa and 1 samples")
})
