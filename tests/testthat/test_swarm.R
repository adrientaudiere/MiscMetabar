suppressWarnings(swarm_error_or_not <-
                   try(system("swarm -h 2>&1", intern = TRUE), silent = TRUE))

if (class(swarm_error_or_not) == "try-error") {
  message(
    "swarm_clustering() and asv2otu(..., method=swarm) can't be
    tested when swarm is not installed"
  )
} else {

  sequences_ex <- c(
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTAATAACGAATTCATTGAATCA",
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTT",
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAGAGGTG",
    "TACCTATGTTGCCTTGGCGGCTAAACCTACC",
    "CGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCATGAGGTG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG"
  )

  test_that("swarm_clustering works fine with phyloseq object", {
    expect_s4_class(data_fungi_swarm <- swarm_clustering(data_fungi), "phyloseq")
    expect_equal(ntaxa(data_fungi_swarm), 1301)
  })

  test_that("swarm_clustering works fine with dna sequences vector", {
    expect_s3_class(sequences_ex_swarm <- swarm_clustering(dna_seq = sequences_ex), "data.frame")
    expect_equal(dim(sequences_ex_swarm), c(12, 10))
  })

  test_that("asv2otu works fine with swarm method", {
    expect_s4_class(
      d_swarm <-
        asv2otu(data_fungi_sp_known, method = "swarm"),
      "phyloseq"
    )

    expect_equal(ntaxa(d_swarm), 600)
    expect_true(nsamples(d_swarm) == nsamples(data_fungi_sp_known))

     expect_s3_class(
       seq_swarm <-
      asv2otu(dna_seq = sequences_ex, method = "swarm"),
      "data.frame"
    )
     expect_equal(dim(seq_swarm), c(12, 10))
  })
}
