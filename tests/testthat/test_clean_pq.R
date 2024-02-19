data_fungi_test <- data_fungi
data_fungi_test@otu_table[, 1] <- rep(0, nrow(data_fungi_test@otu_table))
data_fungi_test@otu_table[10, ] <- rep(0, ncol(data_fungi_test@otu_table))

test_that("clean_pq return a phyloseq object after cleaning empty taxa and samples", {
  expect_s4_class(data_fungi_test, "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test), "phyloseq")
  skip_on_cran()
  expect_s4_class(clean_pq(data_fungi_test, verbose = TRUE), "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test, reorder_asv = TRUE), "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test, rename_asv = TRUE), "phyloseq")
  expect_error(clean_pq(data_fungi_test, force_taxa_as_columns = TRUE, force_taxa_as_rows = TRUE))
})

test_that("clean_pq clean empty taxa and samples", {
  expect_message(clean_pq(data_fungi_test), "Cleaning suppress 2 taxa and 1 samples")
  skip_on_cran()
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
taxa_names(data_fungi_test2) <- paste0("ASV", seq_along(taxa_names(data_fungi_test2)))
test_that("clean_pq works fine with bad taxa_names", {
  expect_s4_class(clean_pq(data_fungi_test2), "phyloseq")
  expect_message(clean_pq(data_fungi_test2), "Cleaning suppress 2 taxa and 1 samples")
})

data_fungi_test3 <- data_fungi_test
taxa_names(data_fungi_test3@refseq) <- paste0("ASV", seq_along(taxa_names(data_fungi_test3)))
test_that("clean_pq works fine with bad taxa_names in refseq", {
  expect_s4_class(clean_pq(data_fungi_test3), "phyloseq")
  expect_message(clean_pq(data_fungi_test3), "Cleaning suppress 2 taxa and 1 samples")
})

data_fungi_test4 <- data_fungi_test
taxa_names(data_fungi_test4@tax_table) <- paste0("ASV", seq_along(taxa_names(data_fungi_test4)))
test_that("clean_pq works fine with bad taxa_names in tax_table", {
  expect_s4_class(clean_pq(data_fungi_test4), "phyloseq")
  expect_message(clean_pq(data_fungi_test4), "Cleaning suppress 2 taxa and 1 samples")
})

data_fungi_test5 <- data_fungi_test
sample_names(data_fungi_test5@tax_table) <- paste0("SAMP", seq_along(taxa_names(data_fungi_test5)))
test_that("clean_pq works fine with bad sample_names in sam_data", {
  expect_s4_class(clean_pq(data_fungi_test5), "phyloseq")
  expect_message(clean_pq(data_fungi_test5), "Cleaning suppress 2 taxa and 1 samples")
})


data_fungi_test6 <- data_fungi_test
sample_names(data_fungi_test6)[1] <- paste0("0", sample_names(data_fungi_test6)[1])
test_that("clean_pq works fine with one sample with a 0 at the start of the name", {
  expect_s4_class(clean_pq(data_fungi_test6), "phyloseq")
  expect_message(expect_message(clean_pq(data_fungi_test6), "Cleaning suppress 2 taxa and 1 samples"))
})
