data(enterotype)

test_that("write_pq function works fine", {
  testFolder <- tempdir()
  expect_silent(write_pq(enterotype, path = testFolder, silent = TRUE))
  expect_message(write_pq(enterotype, path = testFolder))
  expect_message(write_pq(enterotype, one_file = TRUE, path = testFolder))
  expect_s4_class(read_pq(testFolder), "phyloseq")
  new_data_fungi <- read_pq(testFolder)
  expect_equal(ntaxa(new_data_fungi) - ntaxa(enterotype), 0)
  expect_equal(nsamples(new_data_fungi) - nsamples(enterotype), 0)
})
