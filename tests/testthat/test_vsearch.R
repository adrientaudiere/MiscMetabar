skip_on_cran()
sequences_ex <- c(
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG",
  "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCAACCTACCCTGTACACCGCGCCTGGGTCTACCCTCCGGATGACATTTTTAAGACTCTTGTTTTATAGTGAAATTCTGAGTTTTTATACTTAATAAGTTAAAACTTTCAATCTCGGATCTCTTGGCTCTGGCATCGATGAAGAACGCTACGAAATGCTGATAAATAATGTGAATTGCCGAATTCATTGAATCATCGAATCTTTGAACGCACATTGCACCCATTAGTATTCTAGAGTGCATGCCTGTTCCAGCGTCATTTTCAATCCTCAAGCCCCTTATTGCTTGGTGTTGGCAGTTTAGCTGGCTTTATAGTGCTTAACTCCCTAAATATACTGCCTGATTCGCGGTGACCCCAAGCGTAATAATTATTTTCTCGCTTGAGGTG"
)

data("data_fungi", package = "MiscMetabar")
data("data_fungi_sp_known", package = "MiscMetabar")
data("data_fungi_mini", package = "MiscMetabar")
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
    expect_identical(sum(dim(d_vs@otu_table) == dim(d_fast@otu_table)), 2L)
  })

  test_that("vs_search_global works fine with vsearch method", {
    expect_s3_class(
      res <- vs_search_global(
        data_fungi,
        path_to_fasta = "inst/extdata/ex_little.fasta"
      ),
      "data.frame"
    )
    expect_identical(dim(res), c(1420L, 10L))
    expect_s3_class(
      res <-
        vs_search_global(data_fungi, seq2search = sequences_ex),
      "data.frame"
    )
    expect_s3_class(
      res <-
        vs_search_global(
          data_fungi,
          seq2search = Biostrings::DNAStringSet(sequences_ex)
        ),
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

    expect_true(length(chimera_fungi$non_chimera) %in% c(1051, 1088, 1054))
    expect_true(length(chimera_fungi$chimera) %in% c(220, 242, 240))
    expect_true(length(chimera_fungi$borderline) %in% c(112, 127, 126))
  })

  test_that("chimera_detection_vs works fine", {
    expect_s4_class(
      data_fungi_nochim <-
        chimera_removal_vs(data_fungi),
      "phyloseq"
    )

    expect_true(ntaxa(data_fungi_nochim) %in% c(1178, 1200, 1180))
    expect_s4_class(
      data_fungi_nochim_16 <- chimera_removal_vs(
        data_fungi,
        abskew = 16,
        min_seq_length = 10
      ),
      "phyloseq"
    )
    expect_true(ntaxa(data_fungi_nochim_16) %in% c(1259, 1288, 1261))
    expect_s4_class(
      data_fungi_nochim2 <-
        chimera_removal_vs(
          data_fungi,
          type = "Select_only_non_chim_seqlen_filtered"
        ),
      "phyloseq"
    )
    expect_true(ntaxa(data_fungi_nochim2) %in% c(1051, 1088, 1054))

    expect_s4_class(
      data_fungi_chimera <-
        chimera_removal_vs(data_fungi, type = "Select_only_chim"),
      "phyloseq"
    )
    expect_true(ntaxa(data_fungi_chimera) %in% c(220, 242, 240))
  })

  test_that("vsearch_clustering works fine", {
    expect_s4_class(d_vs1 <- vsearch_clustering(data_fungi), "phyloseq")
    expect_identical(ntaxa(d_vs1), 701L)

    expect_s4_class(
      d_vs2 <- vsearch_clustering(
        data_fungi,
        id = 0.98,
        vsearch_cluster_method = "--cluster_size"
      ),
      "phyloseq"
    )
    expect_identical(ntaxa(d_vs2), 817L)

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
    expect_identical(dim(seq_clustered), c(4L, 10L))
  })

  test_that("assign_vsearch_lca works fine", {
    expect_s3_class(
      assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar",
          mustWork = TRUE
        ),
        id = 0.4
      ),
      "tbl_df"
    )

    expect_type(
      assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        id = 0.4,
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
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        lca_cutoff = 0.6
      ),
      "phyloseq"
    )

    data_fungi_mini_new_id90 <- assign_vsearch_lca(
      data_fungi_mini,
      ref_fasta = system.file(
        "extdata",
        "mini_UNITE_fungi.fasta.gz",
        package = "MiscMetabar"
      ),
      behavior = "add_to_phyloseq",
      id = 0.9
    )

    expect_s4_class(data_fungi_mini_new_id90, "phyloseq")

    expect_s4_class(
      assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        behavior = "add_to_phyloseq",
        id = 0.6
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca90 <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        lca_cutoff = 0.9,
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca90_tophit_maxa4 <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
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
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        lca_cutoff = 0.9,
        top_hits_only = TRUE,
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca100_tophit <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        top_hits_only = TRUE,
        behavior = "add_to_phyloseq"
      ),
      "phyloseq"
    )

    expect_s4_class(
      data_fungi_mini_new_lca100_tophit_maxr4 <- assign_vsearch_lca(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        top_hits_only = TRUE,
        behavior = "add_to_phyloseq",
        maxrejects = 4,
        verbose = FALSE
      ),
      "phyloseq"
    )
  })

  test_that("assign_sintax works fine", {
    expect_type(
      assign_sintax(
        data_fungi_mini,
        ref_fasta = system.file(
          "extdata",
          "mini_UNITE_fungi.fasta.gz",
          package = "MiscMetabar"
        ),
        behavior = "return_cmd"
      ),
      "character"
    )

    expect_s4_class(
      data_fungi_mini_new_2 <-
        assign_sintax(
          data_fungi_mini,
          ref_fasta = system.file(
            "extdata",
            "mini_UNITE_fungi.fasta.gz",
            package = "MiscMetabar"
          ),
          min_bootstrap = 0.8,
          behavior = "add_to_phyloseq",
          verbose = FALSE
        ),
      "phyloseq"
    )

    expect_length(
      assignation_results <-
        assign_sintax(
          data_fungi_mini,
          ref_fasta = system.file(
            "extdata",
            "mini_UNITE_fungi.fasta.gz",
            package = "MiscMetabar"
          )
        ),
      2
    )

    # return_matrix must apply min_bootstrap to taxo_value (previously only
    # add_to_phyloseq applied the filter; return_matrix returned raw values)
    res_filt <- assign_sintax(
      data_fungi_mini,
      ref_fasta = system.file(
        "extdata",
        "mini_UNITE_fungi.fasta.gz",
        package = "MiscMetabar"
      ),
      behavior = "return_matrix",
      min_bootstrap = 0.8,
      verbose = FALSE
    )
    rank_cols <- setdiff(names(res_filt$taxo_value), "taxa_names")
    tax_val <- as.matrix(res_filt$taxo_value[, rank_cols, drop = FALSE])
    tax_boot <- as.matrix(res_filt$taxo_bootstrap[, rank_cols, drop = FALSE])
    low_boot <- !is.na(tax_boot) & tax_boot < 0.8
    expect_true(any(low_boot))
    expect_true(all(is.na(tax_val[low_boot])))
  })

  test_that("assign_sintax works with seq2search input", {
    ref <- system.file(
      "extdata",
      "mini_UNITE_fungi.fasta.gz",
      package = "MiscMetabar"
    )
    seqs_to_assign <- refseq(data_fungi_mini)

    # return_cmd works with seq2search and returns a character string
    expect_type(
      assign_sintax(
        seq2search = seqs_to_assign,
        ref_fasta = ref,
        behavior = "return_cmd"
      ),
      "character"
    )

    # return_matrix returns a length-2 list with the expected rank columns
    res <- assign_sintax(
      seq2search = seqs_to_assign,
      ref_fasta = ref,
      verbose = FALSE
    )
    expect_length(res, 2)
    expect_true(all(c("taxo_value", "taxo_bootstrap") %in% names(res)))
    rank_cols <- c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )
    expect_true(all(rank_cols %in% names(res$taxo_value)))
    # one row per input sequence
    expect_equal(nrow(res$taxo_value), length(seqs_to_assign))

    # add_to_phyloseq is not allowed with seq2search
    expect_error(
      assign_sintax(
        seq2search = seqs_to_assign,
        ref_fasta = ref,
        behavior = "add_to_phyloseq"
      ),
      "add_to_phyloseq"
    )
  })

  test_that("assign_sintax return_taxtab returns a ready-to-use matrix", {
    ref <- system.file(
      "extdata",
      "mini_UNITE_fungi.fasta.gz",
      package = "MiscMetabar"
    )
    seqs_to_assign <- refseq(data_fungi_mini)

    tax_mat <- assign_sintax(
      seq2search = seqs_to_assign,
      ref_fasta = ref,
      behavior = "return_taxtab",
      min_bootstrap = 0.8,
      verbose = FALSE
    )
    expect_type(tax_mat, "character")
    expect_true(is.matrix(tax_mat))
    expect_equal(rownames(tax_mat), names(seqs_to_assign))
    rank_cols <- c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )
    expect_equal(colnames(tax_mat), rank_cols)
    expect_equal(nrow(tax_mat), length(seqs_to_assign))

    # min_bootstrap filter is applied (NAs present where bootstrap < 0.8)
    expect_true(any(is.na(tax_mat)))
  })

  test_that("assign_sintax accepts a matrix as seq2search", {
    ref <- system.file(
      "extdata",
      "mini_UNITE_fungi.fasta.gz",
      package = "MiscMetabar"
    )
    seqs <- refseq(data_fungi_mini)
    # Build a dada2-style sequence table: colnames are the DNA sequences
    seqtab <- matrix(1, nrow = 1, ncol = length(seqs))
    colnames(seqtab) <- unname(as.character(seqs))

    tax_mat <- assign_sintax(
      seq2search = seqtab,
      ref_fasta = ref,
      behavior = "return_taxtab",
      verbose = FALSE
    )
    expect_true(is.matrix(tax_mat))
    # rownames are the colnames of the matrix (the sequences)
    expect_equal(rownames(tax_mat), colnames(seqtab))
    expect_equal(nrow(tax_mat), ncol(seqtab))
    rank_cols <- c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    )
    expect_equal(colnames(tax_mat), rank_cols)
  })

  test_that("assign_vsearch_lca return_taxtab with seq2search", {
    ref <- system.file(
      "extdata",
      "mini_UNITE_fungi.fasta.gz",
      package = "MiscMetabar"
    )
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
    tax_mat <- assign_vsearch_lca(
      seq2search = seqs,
      ref_fasta = ref,
      behavior = "return_taxtab",
      verbose = FALSE
    )
    expect_true(is.matrix(tax_mat))
    expect_type(tax_mat, "character")
    expect_equal(colnames(tax_mat), rank_cols)
    expect_equal(rownames(tax_mat), names(seqs))
    expect_equal(nrow(tax_mat), length(seqs))

    # Matrix input (colnames are the sequences)
    tax_mat2 <- assign_vsearch_lca(
      seq2search = seqtab,
      ref_fasta = ref,
      behavior = "return_taxtab",
      verbose = FALSE
    )
    expect_true(is.matrix(tax_mat2))
    expect_equal(colnames(tax_mat2), rank_cols)
    expect_equal(rownames(tax_mat2), colnames(seqtab))
    expect_equal(nrow(tax_mat2), ncol(seqtab))
  })
}
