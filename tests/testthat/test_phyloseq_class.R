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
    expect_s4_class(suppressWarnings(suppressMessages(lulu_pq(data_fungi)))$new_physeq, "phyloseq")
    expect_error(suppressWarnings(suppressMessages(lulu_pq(enterotype)))$new_physeq)
    expect_s4_class(suppressWarnings(suppressMessages(lulu_pq(data_fungi_sp_known, clean_pq = TRUE, verbose = TRUE)))$new_physeq, "phyloseq")
    expect_error(suppressWarnings(suppressMessages(lulu_pq(data_fungi_sp_known, minimum_ratio_type = "avg")))$new_physeq)
  })
}


test_that("lulu works", {
  skip_on_cran()
  suppressWarnings(suppressMessages(res1 <- lulu_pq(data_fungi_sp_known)))
  expect_equal(length(res1), 4)
  expect_equal(ntaxa(res1$new_physeq), 549)
  expect_equal(nsamples(res1$new_physeq), 185)

  suppressWarnings(suppressMessages(res2 <- lulu_pq(data_fungi_sp_known, verbose = TRUE, clean_pq = TRUE)))
  expect_equal(length(res2), 4)
  expect_equal(ntaxa(res2$new_physeq), 549)
  expect_equal(nsamples(res2$new_physeq), 184)
})

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

  data_fungi2 <- data_fungi
  taxa_names(data_fungi2@otu_table) <- paste0("New_names_",taxa_names(data_fungi) )
  expect_error(verify_pq(data_fungi2), "Inconsistency of taxa_names between otu_table and tax_table slots.")

  data_fungi3 <- data_fungi
  sample_names(data_fungi3@sam_data) <- paste0("New_names_", sample_names(data_fungi) )
  expect_error(verify_pq(data_fungi3), "Inconsistency of sample_names between otu_table and sam_data slots.")
})

test_that("verify_tax_table works fine", {
  # Test with clean data and verbose = FALSE (default) - should return silently
  expect_silent(verify_tax_table(data_fungi_mini))

  # Test with clean data and verbose = TRUE - may produce messages/warnings
  expect_silent(suppressMessages(suppressWarnings(verify_tax_table(data_fungi_mini, verbose = TRUE))))

  # Test with NA-like patterns (verbose = TRUE required for warnings)
  data_unclassified <- data_fungi_mini
  data_unclassified@tax_table[1, "Genus"] <- "unclassified_Fungi"
  expect_silent(verify_tax_table(data_unclassified)) # verbose = FALSE
  expect_warning(verify_tax_table(data_unclassified, verbose = TRUE), "NA-like patterns")

  # Test with empty QIIME-style rank
  data_empty_rank <- data_fungi_mini
  data_empty_rank@tax_table[1, "Genus"] <- "g__"
  expect_warning(verify_tax_table(data_empty_rank, verbose = TRUE), "NA-like patterns")

  # Test with short values
  data_short <- data_fungi_mini
  data_short@tax_table[1, "Genus"] <- "sp"
  expect_warning(verify_tax_table(data_short, verbose = TRUE), "less than 4 characters")

  # Test with leading/trailing spaces
  data_spaces <- data_fungi_mini
  data_spaces@tax_table[1, "Genus"] <- " Peziza "
  expect_warning(verify_tax_table(data_spaces, verbose = TRUE), "leading or trailing whitespace")

  # Test with ranks containing only NA
  data_all_na <- data_fungi_mini
  data_all_na@tax_table[, "Species"] <- NA
  expect_warning(verify_tax_table(data_all_na, verbose = TRUE), "only NA values")

  # Test check_taxonomy parameter in verify_pq (implicitly sets verbose = TRUE)
  expect_warning(verify_pq(data_unclassified, check_taxonomy = TRUE), "NA-like patterns")
})

data_fungi_with__P <- data_fungi_mini
data_fungi_with__P@tax_table[, "Phylum"] <- paste0("P__", data_fungi_with__P@tax_table[, "Phylum"])
test_that("simplify_taxo works fine", {
  expect_equal(sum(data_fungi_with__P@tax_table[, "Phylum"] == data_fungi_mini@tax_table[, "Phylum"]), 0)
  expect_equal(simplify_taxo(data_fungi_with__P)@tax_table[, "Phylum"], data_fungi_mini@tax_table[, "Phylum"])
})


test_that("get_file_extension works fine", {
  expect_warning(get_file_extension("inst/extdata/ex_R1_001.fastq.gz"), "There is more than one '.' inside your file path: inst/extdata/ex_R1_001.fastq.gz")
  expect_equal(suppressWarnings(get_file_extension("inst/extdata/ex_R1_001.fastq.gz")), c("fastq", "gz"))
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
  expect_equal(suppressWarnings(count_seq(folder_path = "inst/extdata", pattern = "*.fasta")), c(1000, 500, 3, 2, 5000, 500))
  expect_equal(suppressWarnings(count_seq("inst/extdata/ex_R1_001.fastq.gz")), 2500)
  expect_equal(suppressWarnings(count_seq("inst/extdata/ex_R2_001.fastq.gz")), 2500)
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
