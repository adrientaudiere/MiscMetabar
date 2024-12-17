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

if (!MiscMetabar:::is_vsearch_installed()) {
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
      res <- vs_search_global(data_fungi, path_to_fasta = "inst/extdata/ex_little.fasta"),
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

    expect_true(length(chimera_fungi$non_chimera) %in% c(1051, 1088))
    expect_true(length(chimera_fungi$chimera) %in% c(220, 242))
    expect_true(length(chimera_fungi$borderline) %in% c(112, 127))
  })

  test_that("chimera_detection_vs works fine", {
    expect_s4_class(
      data_fungi_nochim <-
        chimera_removal_vs(data_fungi),
      "phyloseq"
    )

    expect_true(ntaxa(data_fungi_nochim) %in% c(1178, 1200))
    expect_s4_class(
      data_fungi_nochim_16 <- chimera_removal_vs(data_fungi, abskew = 16, min_seq_length = 10),
      "phyloseq"
    )
    expect_true(ntaxa(data_fungi_nochim_16) %in% c(1259, 1288))
    expect_s4_class(
      data_fungi_nochim2 <-
        chimera_removal_vs(data_fungi, type = "Select_only_non_chim"),
      "phyloseq"
    )
    expect_true(ntaxa(data_fungi_nochim2) %in% c(1051, 1088))

    expect_s4_class(
      data_fungi_chimera <-
        chimera_removal_vs(data_fungi, type = "Select_only_chim"),
      "phyloseq"
    )
    expect_true(ntaxa(data_fungi_chimera) %in% c(220, 242))
  })


  test_that("vsearch_clustering works fine", {
    expect_s4_class(d_vs1 <- vsearch_clustering(data_fungi), "phyloseq")
    expect_equal(ntaxa(d_vs1), 701)

    expect_s4_class(
      d_vs2 <- vsearch_clustering(
        data_fungi,
        id = 0.98,
        vsearch_cluster_method = "--cluster_size"
      ),
      "phyloseq"
    )
    expect_equal(ntaxa(d_vs2), 817)

    expect_s4_class(
      d_vs3 <- vsearch_clustering(
        data_fungi,
        id = 0.98,
        vsearch_cluster_method = "--cluster_smallmem",
        vsearch_args = "--strand both --usersort"
      ),
      "phyloseq"
    )


    expect_type(
      seq_clustered <- vsearch_clustering(dna_seq = sequences_ex),
      "list"
    )
    expect_equal(dim(seq_clustered), c(4, 10))
  })

  test_that("assign_vsearch_lca works fine", {
    expect_s3_class(
      assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        id = 0.9
      ),
      "tbl_df"
    )

    expect_type(
      assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        id = 0.9,
        behavior = "return_cmd"
      ),
      "character"
    )

    expect_s4_class(
      data_fungi_mini_new_maxa100 <- assign_vsearch_lca(
        data_fungi_mini,
        top_hits_only = FALSE,
        behavior = "add_to_phyloseq",
        maxaccepts = 5,
        maxreject = 0,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        lca_cutoff = 0.8
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_id90 <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        behavior = "add_to_phyloseq",
        id = 0.9
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca90 <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        lca_cutoff = 0.9,
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca90_tophit_maxa4 <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        lca_cutoff = 0.9,
        top_hits_only = TRUE,
        behavior = "add_to_phyloseq",
        maxaccepts = 4
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca90_tophit <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        lca_cutoff = 0.9,
        top_hits_only = TRUE,
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca100_tophit <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        top_hits_only = TRUE,
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca100_tophit_maxr4 <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        top_hits_only = TRUE,
        behavior = "add_to_phyloseq",
        maxrejects = 4,
        verbose=FALSE
      ),
      "phyloseq"
    )

    res <- lapply(list(
      data_fungi_mini_new_maxa100,
      data_fungi_mini_new_id90,
      data_fungi_mini_new_lca90,
      data_fungi_mini_new_lca90_tophit,
      data_fungi_mini_new_lca90_tophit_maxa4,
      data_fungi_mini_new_lca100_tophit_maxr4,
      data_fungi_mini_new_lca100_tophit
    ), function(el) {
      sum(is.na(el@tax_table[, "G"]))
    })

    expect_equal(res, list(35, 43, 0, 0, 0, 0, 0))
  })


  test_that("assign_sintax works fine", {
    expect_type(
      assign_sintax(
        data_fungi_mini,
        ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
        behavior = "return_cmd"
      ),
      "character"
    )

    expect_s4_class(data_fungi_mini_new <- assign_sintax(data_fungi_mini,
      ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
      behavior = "add_to_phyloseq"
    ), "phyloseq")

    expect_s4_class(data_fungi_mini_new_2 <-
      assign_sintax(
        data_fungi_mini,
        ref_fasta = system.file("extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        min_boostrap = 0.8,
        behavior = "add_to_phyloseq",
        verbose=FALSE
      ), "phyloseq")

    expect_length(assignation_results <-
      assign_sintax(data_fungi_mini,
        ref_fasta = system.file("extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        )
      ), 2)
  })
}
