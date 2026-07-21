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
  message("denoised_reads() can't be tested when vsearch is not installed")
} else {
  test_that("denoised_reads works with vsearch UNOISE (dna_seq and files)", {
    pq <- denoised_reads(dna_seq = sequences_ex)
    expect_s4_class(pq, "phyloseq")
    expect_equal(nsamples(pq), 1)
    expect_equal(sum(taxa_sums(pq)), length(sequences_ex))
    expect_true(all(grepl("^ASV_", taxa_names(pq))))
    expect_identical(taxa_names(pq), names(phyloseq::refseq(pq)))

    pq2 <- denoised_reads(path_to_fastx = c(f1, f2))
    expect_s4_class(pq2, "phyloseq")
    expect_equal(nsamples(pq2), 2)
    expect_setequal(sample_names(pq2), c("ex_R1_001", "ex_R1_002"))
    expect_true(all(grepl("^ASV_", taxa_names(pq2))))
    expect_identical(
      as.character(sample_data(pq2)$source_file),
      c(f1, f2)
    )
  })

  test_that("denoised_reads honours minsize and ;size= annotations", {
    sized <- sequences_ex[!duplicated(sequences_ex)]
    names(sized) <- paste0(
      "s",
      seq_along(sized),
      ";size=",
      c(5, 4, 3, 2, 1, 1, 3)
    )

    pq <- denoised_reads(dna_seq = sized)
    expect_equal(sum(taxa_sums(pq)), 19)

    pq_min <- denoised_reads(dna_seq = sized, minsize = 2)
    expect_equal(sum(taxa_sums(pq_min)), 17)

    expect_error(
      denoised_reads(dna_seq = sized, minsize = 100),
      "No sequence left"
    )
  })

  test_that("denoised_reads validates its inputs", {
    expect_error(
      denoised_reads(path_to_fastx = f1, dna_seq = sequences_ex),
      "exactly one"
    )
    expect_error(denoised_reads(), "exactly one")
    expect_error(
      denoised_reads(path_to_fastx = file.path(tempdir(), "nope.fasta")),
      "not found"
    )
    expect_error(
      denoised_reads(path_to_fastx = c(f1, f1)),
      "not unique"
    )
  })

  test_that("denoised_reads with dada2 requires fastq input", {
    expect_error(
      denoised_reads(dna_seq = sequences_ex, method = "dada2"),
      "requires fastq"
    )
    expect_error(
      denoised_reads(path_to_fastx = f1, method = "dada2"),
      "requires fastq"
    )
  })
}

if (
  MiscMetabar:::is_vsearch_installed() && MiscMetabar:::is_swarm_installed()
) {
  test_that("denoised_reads works with swarm d = 1", {
    pq <- denoised_reads(dna_seq = sequences_ex, method = "swarm")
    expect_s4_class(pq, "phyloseq")
    expect_equal(sum(taxa_sums(pq)), length(sequences_ex))
    expect_true(all(grepl("^ASV_", taxa_names(pq))))
  })

  test_that("denoised_reads warns when swarm d > 1", {
    expect_warning(
      denoised_reads(dna_seq = sequences_ex, method = "swarm", d = 3),
      "cluster_reads"
    )
  })
} else {
  message("swarm method of denoised_reads() requires the swarm binary")
}

test_that("denoised_reads works with the classical dada2 pipeline", {
  skip_if_not_installed("dada2")
  fq1 <- system.file("extdata", "ex_R1_001.fastq.gz", package = "MiscMetabar")
  fq2 <- system.file("extdata", "ex.fastq", package = "MiscMetabar")

  suppressMessages(
    pq <- denoised_reads(path_to_fastx = c(fq1, fq2), method = "dada2")
  )
  expect_s4_class(pq, "phyloseq")
  expect_equal(nsamples(pq), 2)
  expect_setequal(sample_names(pq), c("ex_R1_001", "ex"))
  expect_true(all(grepl("^ASV_", taxa_names(pq))))
  expect_identical(taxa_names(pq), names(phyloseq::refseq(pq)))
  expect_equal(ntaxa(pq), length(phyloseq::refseq(pq)))

  suppressMessages(
    pq_nochim <- denoised_reads(
      path_to_fastx = fq1,
      method = "dada2",
      remove_chimeras = FALSE
    )
  )
  expect_s4_class(pq_nochim, "phyloseq")
  expect_equal(nsamples(pq_nochim), 1)
})
