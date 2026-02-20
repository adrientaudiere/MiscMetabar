data(data_fungi)

test_that("assign_blastn works", {
  skip_on_cran()
  ref_fasta <- Biostrings::readDNAStringSet(system.file(
    "extdata",
    "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar",
    mustWork = TRUE
  ))

  # assign_blastn(data_fungi_mini, ref_fasta = ref_fasta) # error because not
  # enough sequences in db so none blast query passed the filters.
  # So we used low score filter hereafter.

  mat <- assign_blastn(
    data_fungi_mini,
    ref_fasta = ref_fasta,
    method_algo = "top-hit",
    min_id = 70,
    min_e_value = 1e-3,
    min_cover = 50,
    min_bit_score = 20
  )

  expect_equal(dim(mat), c(41, 9))

  expect_error(assign_blastn(dna, database = "non_existent_db"))
})


data_fungi_mini_2asv <- subset_taxa_pq(
  data_fungi_mini,
  taxa_names(data_fungi_mini) %in% c("ASV7", "ASV8"),
  taxa_names_from_physeq = TRUE
)

test_that("assign_dada2 works", {
  if (requireNamespace("stringr", quietly = TRUE)) {
    data_fungi_mini2 <- assign_dada2(
      data_fungi_mini_2asv,
      ref_fasta = system.file(
        "extdata",
        "mini_UNITE_fungi.fasta.gz",
        package = "MiscMetabar"
      ),
      suffix = "_dada2",
      from_sintax = TRUE
    )
    expect_s4_class(data_fungi_mini2, "phyloseq")
    expect_equal(ncol(data_fungi_mini2@tax_table), 20)
  }
})

test_that("assign_idtaxa works", {
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
    data_fungi_mini2 <- assign_idtaxa(
      data_fungi_mini_2asv,
      fasta_for_training = system.file(
        "extdata",
        "mini_UNITE_fungi.fasta.gz",
        package = "MiscMetabar"
      ),
      threshold = 20,
      behavior = "add_to_phyloseq"
    )
  }
  expect_s4_class(data_fungi_mini2, "phyloseq")
  expect_equal(ncol(data_fungi_mini2@tax_table), 16)
})
