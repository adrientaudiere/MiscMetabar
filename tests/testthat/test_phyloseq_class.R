data(enterotype)
data("data_fungi_sp_known")

sequences_ex <- c(
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG"
)
test_that("asv2otu works fine with Clusterize method", {
  expect_s4_class(data_fungi_sp_known, "phyloseq")
  expect_s4_class(asv2otu(data_fungi_sp_known), "phyloseq")

  expect_error(asv2otu(enterotype, seq_names = sequences_ex))
  expect_error(asv2otu(enterotype, method = "vsearch"))
  expect_error(asv2otu(enterotype))
  expect_error(asv2otu())
  expect_error(asv2otu(method = "VsaerCh"))
  expect_error(asv2otu(data_fungi_sp_known, seq_names = sequences_ex))
})




suppressWarnings(vsearch_error_or_not <- try(system("vsearch 2>&1", intern = TRUE), silent = TRUE))

if (class(vsearch_error_or_not) == "try-error") {
  message("lulu_phyloseq() can't be tested when vsearch is not installed")
} else {
  test_that("lulu_pq works fine", {
    expect_s4_class(lulu_pq(data_fungi_sp_known)$new_physeq, "phyloseq")
    expect_error(lulu_pq(enterotype)$new_physeq, "phyloseq")
    expect_s4_class(lulu_pq(data_fungi_sp_known, clean_pq = TRUE, verbose = TRUE)$new_physeq, "phyloseq")
  })
}


test_that("as_binary_otu_table works fine", {
  expect_s4_class(as_binary_otu_table(data_fungi_sp_known), "phyloseq")
  expect_s4_class(as_binary_otu_table(enterotype), "phyloseq")
})

data_fungi_taxaSeq <- data_fungi
taxa_names(data_fungi_taxaSeq) <- as.character(data_fungi_taxaSeq@refseq)
data_fungi_taxaSeq@refseq <- NULL
test_that("add_dna_to_phyloseq works fine", {
  expect_silent(data_fungi_taxaSeq <- add_dna_to_phyloseq(data_fungi_taxaSeq))
  expect_equal(taxa_names(data_fungi_taxaSeq), paste0("ASV_", seq_len(length(taxa_names(data_fungi_taxaSeq)))))
})


df <- data_fungi
df <- unclass(data_fungi)
test_that("verify_pq works fine", {
  expect_error(verify_pq(df), "The physeq argument is not a valid phyloseq object.")
})

data_fungi_with__P <- data_fungi
data_fungi_with__P@tax_table[, "Phylum"] <- paste0("P__", data_fungi_with__P@tax_table[, "Phylum"])
test_that("simplify_taxo works fine", {
  expect_equal(sum(data_fungi_with__P@tax_table[, "Phylum"] == data_fungi@tax_table[, "Phylum"]), 0)
  expect_equal(simplify_taxo(data_fungi_with__P)@tax_table[, "Phylum"], data_fungi@tax_table[, "Phylum"])
})


test_that("get_file_extension works fine", {
  expect_equal(get_file_extension("inst/extdata/ex_R1_001.fastq.gz"), c("fastq", "gz"))
  expect_equal(get_file_extension("inst/extdata/ex.fasta"), "fasta")
})


test_that("perc works fine", {
  expect_equal(perc(20, 200), 10)
  expect_equal(perc(20, 200, add_symbol = TRUE), "10%")
  expect_equal(perc(0.1), 10)
  expect_equal(perc(0.1, add_symbol = TRUE), "10%")
})

test_that("count_seq works fine", {
  expect_equal(count_seq("inst/extdata/ex_R1_001.fastq.gz"), 2500)
  expect_equal(count_seq("inst/extdata/ex_R2_001.fastq.gz"), 2500)
  expect_equal(count_seq("inst/extdata/ex.fasta"), 3)
  expect_equal(count_seq("inst/extdata/ex.fastq"), 4)
  expect_error(count_seq("tests/testthat.R"), "The file extension R is not supported.")
})

# test_that("add_new_taxonomy_pq works fine", {
#   expect_s4_class(add_new_taxonomy_pq(data_fungi, "inst/extdata/1000_sp_UNITE_sh_general_release_dynamic.fasta"), "phyloseq")
# })
