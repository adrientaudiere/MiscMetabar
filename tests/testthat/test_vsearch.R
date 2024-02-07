sequences_ex <- c(
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG"
)

data("data_fungi")
df_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
df_basidio <-
  subset_taxa_pq(df_basidio, colSums(df_basidio@otu_table) > 1000)
# path_db <- "inst/extdata/100_sp_UNITE_sh_general_release_dynamic.fasta"

suppressWarnings(vsearch_error_or_not <-
                   try(system("vsearch 2>&1", intern = TRUE), silent = TRUE))

if (class(vsearch_error_or_not) == "try-error") {
  message(
    "vs_search_global() and asv2otu(..., method=vsearch) can't be tested when vsearch is not installed"
  )
} else {
  test_that("asv2otu works fine with vsearch method", {
    expect_s4_class(
      d_vs <-
        asv2otu(data_fungi_sp_known, method = "vsearch"),
      "phyloseq"
    )
    expect_s4_class(
      d_fast <- asv2otu(
        data_fungi_sp_known,
        method = "vsearch",
        vsearch_cluster_method = "--cluster_fast"
      ),
      "phyloseq"
    )
    expect_s3_class(
      asv2otu(dna_seq = sequences_ex, method = "vsearch"),
      "data.frame"
    )
    expect_true(sum(!d_fast@refseq == d_vs@refseq) > 0)
    expect_equal(sum(dim(d_vs@otu_table) == dim(d_fast@otu_table)), 2)
  })

  test_that("vs_search_global works fine with vsearch method", {
    expect_s3_class(
      res <- vs_search_global(data_fungi,
        path_to_fasta = "inst/extdata/ex_little.fasta"
      ),
      "data.frame"
    )
    expect_equal(dim(res), c(1420, 10))
    expect_s3_class(
      res <-
        vs_search_global(data_fungi, sequences_ex),
      "data.frame"
    )
    expect_s3_class(
      res <-
        vs_search_global(data_fungi, Biostrings::DNAStringSet(sequences_ex)),
      "data.frame"
    )
  })

  test_that("chimera_detection_vs works fine", {
    expect_type(
      chimera_fungi <- chimera_detection_vs(
        seq2search = data_fungi@refseq,
        nb_seq = taxa_sums(data_fungi)
      ),
      "list"
    )
    expect_s4_class(chimera_fungi$non_chimera, "AAStringSet")

    expect_equal(length(chimera_fungi$non_chimera), 1051)
    expect_equal(length(chimera_fungi$chimera), 242)
    expect_equal(length(chimera_fungi$borderline), 127)
  })

  test_that("chimera_detection_vs works fine", {
    expect_s4_class(
      data_fungi_nochim <-
        chimera_removal_vs(data_fungi),
      "phyloseq"
    )
    expect_equal(ntaxa(data_fungi_nochim), 1178)
    expect_s4_class(
      data_fungi_nochim_16 <- chimera_removal_vs(data_fungi,
        abskew = 16,
        min_seq_length = 10
      ),
      "phyloseq"
    )
    expect_equal(ntaxa(data_fungi_nochim_16), 1259)
    expect_s4_class(
      data_fungi_nochim2 <-
        chimera_removal_vs(data_fungi, type = "Select_only_non_chim"),
      "phyloseq"
    )
    expect_equal(ntaxa(data_fungi_nochim2), 1051)

    expect_s4_class(
      data_fungi_chimera <-
        chimera_removal_vs(data_fungi, type = "Select_only_chim"),
      "phyloseq"
    )
    expect_equal(ntaxa(data_fungi_chimera), 242)
  })


  test_that("vsearch_clustering works fine", {
    expect_s4_class(d_vs1 <- vsearch_clustering(data_fungi), "phyloseq")
    expect_equal(ntaxa(d_vs1), 701)

    expect_s4_class(d_vs2 <- vsearch_clustering(data_fungi,
      id = 0.98,
      vsearch_cluster_method = "--cluster_size"
    ), "phyloseq")
    expect_equal(ntaxa(d_vs2), 817)

    expect_s4_class(d_vs3 <- vsearch_clustering(data_fungi,
                                               id = 0.98,
                                               vsearch_cluster_method = "--cluster_smallmem",
                                               vsearch_args = "--strand both --usersort"
    ), "phyloseq")


    expect_type(seq_clustered <- vsearch_clustering(dna_seq = sequences_ex),
                "list")
    expect_equal(dim(seq_clustered), c(4, 10))
  })
}
