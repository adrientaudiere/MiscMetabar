data(enterotype)

sequences_ex <- c(
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG"
)
test_that("asv2otu works fine with Clusterize method", {
  expect_s4_class(data_fungi_mini, "phyloseq")
  skip_on_cran()
  expect_s4_class(suppressWarnings(asv2otu(data_fungi_mini)), "phyloseq")
  expect_s4_class(suppressWarnings(asv2otu(data_fungi_mini,
    method_clusterize = "longest"
  )), "phyloseq")

  expect_error(asv2otu(enterotype, dna_seq = sequences_ex))
  expect_error(asv2otu(enterotype, method = "vsearch"))
  expect_error(asv2otu(enterotype))
  expect_error(asv2otu())
  expect_error(asv2otu(data_fungi_mini, method = "VsaerCh"))
  expect_error(asv2otu(data_fungi_mini, dna_seq = sequences_ex))
})

if (!is_vsearch_installed()) {
  message("lulu_pq() can't be tested when vsearch is not installed")
} else {
  test_that("lulu_pq works fine", {
    skip_on_cran()
    expect_s4_class(lulu_pq(data_fungi)$new_physeq, "phyloseq")
    expect_error(lulu_pq(enterotype)$new_physeq)
    expect_s4_class(lulu_pq(data_fungi_sp_known, clean_pq = TRUE, verbose = TRUE)$new_physeq, "phyloseq")
    expect_error(lulu_pq(data_fungi_sp_known, minimum_ratio_type = "avg")$new_physeq)
  })
}

suppressWarnings(mumu_error_or_not <- try(system("mumu --help", intern = TRUE), silent = TRUE))

if (class(mumu_error_or_not) == "try-error") {
  message("mumu_pq() can't be tested when mumu is not installed")
} else {
  test_that("mumu_pq works fine", {
    skip_on_cran()
    expect_s4_class(mumu_pq(data_fungi_mini)$new_physeq, "phyloseq")
    expect_error(mumu_pq(enterotype)$new_physeq, "phyloseq")
    expect_s4_class(mumu_pq(data_fungi_mini, clean_pq = TRUE, verbose = TRUE)$new_physeq, "phyloseq")
  })
}


test_that("as_binary_otu_table works fine", {
  skip_on_cran()
  expect_s4_class(as_binary_otu_table(data_fungi_mini), "phyloseq")
  expect_s4_class(as_binary_otu_table(enterotype), "phyloseq")
})

data_fungi_taxaSeq <- data_fungi_mini
taxa_names(data_fungi_taxaSeq) <- as.character(data_fungi_taxaSeq@refseq)
data_fungi_taxaSeq@refseq <- NULL
test_that("add_dna_to_phyloseq works fine", {
  expect_silent(data_fungi_taxaSeq <- add_dna_to_phyloseq(data_fungi_taxaSeq))
  expect_equal(taxa_names(data_fungi_taxaSeq), paste0("Taxa_", seq_len(length(taxa_names(data_fungi_taxaSeq)))))
})


test_that("verify_pq works fine", {
  expect_error(verify_pq(unclass(data_fungi_mini)), "The physeq argument is not a valid phyloseq object.")
  expect_silent(suppressWarnings(verify_pq(data_fungi_mini, verbose = TRUE)))
})

data_fungi_with__P <- data_fungi_mini
data_fungi_with__P@tax_table[, "Phylum"] <- paste0("P__", data_fungi_with__P@tax_table[, "Phylum"])
test_that("simplify_taxo works fine", {
  expect_equal(sum(data_fungi_with__P@tax_table[, "Phylum"] == data_fungi_mini@tax_table[, "Phylum"]), 0)
  expect_equal(simplify_taxo(data_fungi_with__P)@tax_table[, "Phylum"], data_fungi_mini@tax_table[, "Phylum"])
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
  skip_on_os("windows")
  skip_on_os("mac")
  skip_on_cran()
  expect_equal(count_seq(folder_path = "inst/extdata", pattern = "*.fasta"), c(1000, 3, 2))
  expect_equal(count_seq("inst/extdata/ex_R1_001.fastq.gz"), 2500)
  expect_equal(count_seq("inst/extdata/ex_R2_001.fastq.gz"), 2500)
  expect_equal(count_seq("inst/extdata/ex.fasta"), 3)
  expect_equal(count_seq("inst/extdata/ex.fastq"), 4)
  expect_error(count_seq("tests/testthat.R"), "The file extension R is not supported.")
})

# test_that("add_new_taxonomy_pq works fine", {
#   expect_s4_class(add_new_taxonomy_pq(data_fungi, "inst/extdata/100_sp_UNITE_sh_general_release_dynamic.fasta"), "phyloseq")
# })

test_that("reorder_taxa_pq works fine", {
  expect_silent(data_fungi_asc_ordered_by_abundance <- reorder_taxa_pq(
    data_fungi,
    taxa_names(data_fungi)[order(taxa_sums(data_fungi))]
  ))
  expect_equal(unclass(data_fungi_asc_ordered_by_abundance@tax_table[, "Genus"])[2], "Peziza")
})
