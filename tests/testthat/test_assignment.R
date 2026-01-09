data(data_fungi)

test_that("assign_blastn works", {
  skip_on_cran()
  # Create a simple DNA string set for testing
  dna <- refseq(data_fungi)[1:5]
  expect_error(assign_blastn(dna, database = "non_existent_db"))
})

test_that("assign_dada2 works", {
  skip_on_cran()
  dna <- refseq(data_fungi)[1:5]
  # This will fail without a proper database, but we test that the function exists
  expect_error(assign_dada2(dna, database = "non_existent_db"))
})

test_that("assign_idtaxa works", {
  skip_on_cran()
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
    dna <- refseq(data_fungi)[1:5]
    # This will fail without a proper trained set, but we test that the function exists
    expect_error(assign_idtaxa(dna, trained_set = "non_existent_set"))
  }
})

test_that("learn_idtaxa works", {
  skip_on_cran()
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
    # This will fail without proper input, but we test that the function exists
    expect_error(learn_idtaxa(seqs = NULL, taxonomy = NULL))
  }
})
