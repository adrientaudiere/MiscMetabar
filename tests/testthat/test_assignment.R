skip_on_cran()
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

  expect_identical(dim(mat), c(41L, 9L))

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
    expect_identical(ncol(data_fungi_mini2@tax_table), 20L)
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
  expect_identical(ncol(data_fungi_mini2@tax_table), 16L)
})

test_that("assign_blastn return_taxtab with seq2search", {
  ref_fasta <- Biostrings::readDNAStringSet(system.file(
    "extdata",
    "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar",
    mustWork = TRUE
  ))
  seqs <- refseq(data_fungi_mini)
  seqtab <- matrix(1, nrow = 1, ncol = length(seqs))
  colnames(seqtab) <- unname(as.character(seqs))
  rank_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )

  # DNAStringSet input
  tax_mat <- assign_blastn(
    seq2search = seqs,
    ref_fasta = ref_fasta,
    method_algo = "top-hit",
    min_id = 70,
    min_e_value = 1e-3,
    min_cover = 50,
    min_bit_score = 20,
    behavior = "return_taxtab"
  )
  expect_true(is.matrix(tax_mat))
  expect_type(tax_mat, "character")
  expect_equal(colnames(tax_mat), rank_cols)
  expect_true(nrow(tax_mat) > 0)

  # Matrix input (colnames are the sequences)
  tax_mat2 <- assign_blastn(
    seq2search = seqtab,
    ref_fasta = ref_fasta,
    method_algo = "top-hit",
    min_id = 70,
    min_e_value = 1e-3,
    min_cover = 50,
    min_bit_score = 20,
    behavior = "return_taxtab"
  )
  expect_true(is.matrix(tax_mat2))
  expect_equal(colnames(tax_mat2), rank_cols)
  expect_true(all(rownames(tax_mat2) %in% colnames(seqtab)))

  # add_to_phyloseq is not allowed with seq2search
  expect_error(
    assign_blastn(
      seq2search = seqs,
      ref_fasta = ref_fasta,
      behavior = "add_to_phyloseq"
    ),
    "add_to_phyloseq"
  )
})

test_that("assign_idtaxa return_taxtab with seq2search", {
  skip_if_not_installed("DECIPHER")
  ref_path <- system.file(
    "extdata",
    "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar"
  )
  seqs <- refseq(data_fungi_mini_2asv)
  seqtab <- matrix(1, nrow = 1, ncol = length(seqs))
  colnames(seqtab) <- unname(as.character(seqs))
  rank_cols <- c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  )

  # DNAStringSet input
  tax_mat <- assign_idtaxa(
    seq2search = seqs,
    fasta_for_training = ref_path,
    threshold = 20,
    behavior = "return_taxtab",
    verbose = FALSE
  )
  expect_true(is.matrix(tax_mat))
  expect_type(tax_mat, "character")
  expect_equal(colnames(tax_mat), rank_cols)
  expect_equal(rownames(tax_mat), names(seqs))
  expect_equal(nrow(tax_mat), length(seqs))

  # Matrix input (colnames are the sequences)
  tax_mat2 <- assign_idtaxa(
    seq2search = seqtab,
    fasta_for_training = ref_path,
    threshold = 20,
    behavior = "return_taxtab",
    verbose = FALSE
  )
  expect_true(is.matrix(tax_mat2))
  expect_equal(colnames(tax_mat2), rank_cols)
  expect_equal(rownames(tax_mat2), colnames(seqtab))

  # add_to_phyloseq is not allowed with seq2search
  expect_error(
    assign_idtaxa(
      seq2search = seqs,
      fasta_for_training = ref_path,
      threshold = 20,
      behavior = "add_to_phyloseq"
    ),
    "add_to_phyloseq"
  )
})
