data("data_fungi")

data_fungi_test <- data_fungi
data_fungi_test@otu_table[ ,1] <- rep(0, nrow(data_fungi_test@otu_table))
data_fungi_test@otu_table[10,] <- rep(0, ncol(data_fungi_test@otu_table))

test_that("clean_pq return a phyloseq object after cleaning empty taxa and samples", {
  expect_s4_class(data_fungi_test, "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test), "phyloseq")
})

test_that("clean_pq clean empty taxa and samples", {
  expect_message(clean_pq(data_fungi_test), "Supress 2 taxa and 1 samples")
  expect_equal(nrow(data_fungi_test@otu_table) - 
      nrow(clean_pq(data_fungi_test)@otu_table), 1)
  expect_equal(ncol(data_fungi_test@otu_table) - 
      ncol(clean_pq(data_fungi_test)@otu_table), 2)
})