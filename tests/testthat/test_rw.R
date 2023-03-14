data(data_fungi)

test_that("write_pq function works fine", {
    testFolder <- tempdir()
    expect_silent(write_pq(data_fungi, path = testFolder))
    expect_snapshot(write_pq(data_fungi, path = testFolder))
    expect_s4_class(read_pq(testFolder), "phyloseq")
    new_data_fungi <- read_pq(testFolder)
    expect_equal(ntaxa(new_data_fungi) - ntaxa(data_fungi), 0)
    expect_equal(nsamples(new_data_fungi) - nsamples(data_fungi), 0)
})

