skip_on_cran()

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

f1 <- system.file("extdata", "ex_R1_001.fasta", package = "MiscMetabar")
f2 <- system.file("extdata", "ex_R1_002.fasta", package = "MiscMetabar")

if (!MiscMetabar:::is_vsearch_installed()) {
  message("cluster_reads() can't be tested when vsearch is not installed")
} else {
  test_that("cluster_reads works with vsearch (dna_seq and files)", {
    pq <- cluster_reads(dna_seq = sequences_ex, id = 0.9)
    expect_s4_class(pq, "phyloseq")
    expect_equal(nsamples(pq), 1)
    expect_equal(sum(taxa_sums(pq)), length(sequences_ex))
    expect_true(all(grepl("^OTU_", taxa_names(pq))))
    expect_identical(taxa_names(pq), names(phyloseq::refseq(pq)))
    expect_false(is.unsorted(rev(rowSums(otu_table(pq)))))

    pq2 <- cluster_reads(path_to_fastx = c(f1, f2))
    expect_s4_class(pq2, "phyloseq")
    expect_equal(nsamples(pq2), 2)
    expect_setequal(sample_names(pq2), c("ex_R1_001", "ex_R1_002"))
    expect_true(all(grepl("^OTU_", taxa_names(pq2))))
    expect_identical(
      as.character(sample_data(pq2)$source_file),
      c(f1, f2)
    )
  })

  test_that("cluster_reads validates its inputs", {
    expect_error(
      cluster_reads(path_to_fastx = f1, dna_seq = sequences_ex),
      "exactly one"
    )
    expect_error(cluster_reads(), "exactly one")
    expect_error(
      cluster_reads(path_to_fastx = file.path(tempdir(), "nope.fasta")),
      "not found"
    )
    expect_error(
      cluster_reads(path_to_fastx = c(f1, f1)),
      "not unique"
    )
  })

  test_that("cluster_reads has no dada2 method", {
    expect_error(
      cluster_reads(dna_seq = sequences_ex, method = "dada2"),
      "should be one of"
    )
  })
}

if (
  MiscMetabar:::is_vsearch_installed() && MiscMetabar:::is_swarm_installed()
) {
  test_that("cluster_reads works with swarm d > 1", {
    pq <- cluster_reads(dna_seq = sequences_ex, method = "swarm", d = 3)
    expect_s4_class(pq, "phyloseq")
    expect_equal(sum(taxa_sums(pq)), length(sequences_ex))
    expect_true(all(grepl("^OTU_", taxa_names(pq))))
  })

  test_that("cluster_reads warns when swarm d = 1", {
    expect_warning(
      cluster_reads(dna_seq = sequences_ex, method = "swarm", d = 1),
      "denoised_reads"
    )
  })
} else {
  message("swarm method of cluster_reads() requires the swarm binary")
}
