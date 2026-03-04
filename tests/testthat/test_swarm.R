if (!MiscMetabar:::is_swarm_installed()) {
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
    expect_s4_class(
      data_fungi_swarm <- swarm_clustering(data_fungi),
      "phyloseq"
    )
    expect_identical(ntaxa(data_fungi_swarm), 1301L)
  })

  test_that("swarm_clustering works fine with dna sequences vector", {
    expect_s3_class(
      sequences_ex_swarm <- swarm_clustering(dna_seq = sequences_ex),
      "data.frame"
    )
    expect_identical(dim(sequences_ex_swarm), c(12L, 10L))
  })

  test_that("asv2otu works fine with swarm method", {
    expect_s4_class(
      d_swarm <-
        asv2otu(data_fungi_sp_known, method = "swarm"),
      "phyloseq"
    )

    expect_identical(ntaxa(d_swarm), 600L)
    expect_true(nsamples(d_swarm) == nsamples(data_fungi_sp_known))

    expect_s3_class(
      seq_swarm <-
        asv2otu(dna_seq = sequences_ex, method = "swarm"),
      "data.frame"
    )
    expect_identical(dim(seq_swarm), c(12L, 10L))
  })
}
