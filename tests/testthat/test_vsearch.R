sequences_ex <- c(
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG"
)

data("data_fungi")
df_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
df_basidio <- subset_taxa_pq(df_basidio, colSums(df_basidio@otu_table) > 1000)
path_db <- "inst/extdata/1000_sp_UNITE_sh_general_release_dynamic.fasta"

suppressWarnings(vsearch_error_or_not <- try(system("vsearch 2>&1", intern = TRUE), silent = TRUE))

if (class(vsearch_error_or_not) == "try-error") {
  message("vs_search_global() and asv2otu(..., method=vsearch) can't be tested when vsearch is not installed")
} else {
  test_that("asv2otu works fine with vsearch method", {
    expect_s4_class(d_vs <- asv2otu(data_fungi_sp_known, method = "vsearch"), "phyloseq")
    expect_s4_class(d_fast <- asv2otu(data_fungi_sp_known, method = "vsearch", vsearch_cluster_method = "--cluster_fast"), "phyloseq")
    expect_s3_class(asv2otu(seq_names = sequences_ex, method = "vsearch"), "data.frame")
    expect_true(sum(!d_fast@refseq == d_vs@refseq) > 0)
    expect_equal(sum(dim(d_vs@otu_table) == dim(d_fast@otu_table)), 2)
  })

  test_that("vs_search_global works fine with vsearch method", {
    expect_s3_class(res <- vs_search_global(data_fungi, path_to_fasta = "inst/extdata/ex_little.fasta"), "data.frame")
    expect_equal(dim(res), c(1420, 10))
    expect_s3_class(res <- vs_search_global(data_fungi, sequences_ex), "data.frame")
    expect_s3_class(res <- vs_search_global(data_fungi, Biostrings::DNAStringSet(sequences_ex)), "data.frame")
  })
}