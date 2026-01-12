data(data_fungi)

test_that("assign_blastn works", {
  skip_on_cran()
  ref_fasta <- Biostrings::readDNAStringSet(system.file("extdata",
  "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar", mustWork = TRUE
))

# assign_blastn(data_fungi_mini, ref_fasta = ref_fasta) # error because not
# enough sequences in db so none blast query passed the filters.
# So we used low score filter hereafter.

mat <- assign_blastn(data_fungi_mini,
  ref_fasta = ref_fasta,
  method_algo = "top-hit", min_id = 70, min_e_value = 1e-3, min_cover = 50,
  min_bit_score = 20
)
  
expect_equal(dim(mat), c(41, 9))
  
  expect_error(assign_blastn(dna, database = "non_existent_db"))
})

test_that("assign_dada2 works", {
  #TODO
})

test_that("assign_idtaxa works", {
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
  #TODO
  }
})

test_that("learn_idtaxa works", {
  if (requireNamespace("DECIPHER", quietly = TRUE)) {
  #TODO
  }
})
