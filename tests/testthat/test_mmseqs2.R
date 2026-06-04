data(data_fungi)

test_that("find_mmseqs2 returns a string", {
  expect_type(find_mmseqs2(), "character")
})

test_that("is_mmseqs2_installed returns logical", {
  expect_type(is_mmseqs2_installed(), "logical")
})

test_that("assign_mmseqs2 validates inputs", {
  skip_on_cran()
  expect_error(
    assign_mmseqs2(
      database = file.path(tempdir(), "fake_db"),
      physeq = data_fungi_mini,
      lca_ranks = c("genus", "species"),
      column_names = c("Genus")
    ),
    "same length"
  )

  expect_error(
    assign_mmseqs2(physeq = data_fungi_mini),
    "ref_fasta.*database"
  )

  expect_error(
    assign_mmseqs2(
      physeq = data_fungi_mini,
      ref_fasta = "some.fasta",
      database = "some_db"
    ),
    "both"
  )
})

test_that("assign_mmseqs2 works with ref_fasta", {
  skip_on_cran()
  skip_if_not(is_mmseqs2_installed(), "MMseqs2 is not installed")

  ref_fasta <- Biostrings::readDNAStringSet(system.file(
    "extdata",
    "mini_UNITE_fungi.fasta.gz",
    package = "MiscMetabar",
    mustWork = TRUE
  ))

  res <- assign_mmseqs2(data_fungi_mini, ref_fasta = ref_fasta)
  expect_true("taxa_names" %in% colnames(res))
  expect_true("Kingdom_mmseqs2" %in% colnames(res))
  expect_true(any(!is.na(res$Kingdom_mmseqs2)))
})

test_that("mmseqs2_clustering works", {
  skip_on_cran()
  skip_if_not(is_mmseqs2_installed(), "MMseqs2 is not installed")

  d_mm <- mmseqs2_clustering(data_fungi_mini)
  expect_s4_class(d_mm, "phyloseq")
  expect_lte(phyloseq::ntaxa(d_mm), phyloseq::ntaxa(data_fungi_mini))
})

test_that("postcluster_pq works with mmseqs2 method", {
  skip_on_cran()
  skip_if_not(is_mmseqs2_installed(), "MMseqs2 is not installed")

  d_mm <- postcluster_pq(data_fungi_mini, method = "mmseqs2")
  expect_s4_class(d_mm, "phyloseq")
  expect_lte(phyloseq::ntaxa(d_mm), phyloseq::ntaxa(data_fungi_mini))
})
