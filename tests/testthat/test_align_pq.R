skip_on_cran()

data("data_fungi_mini")

df <- suppressMessages(subset_taxa_pq(
  data_fungi_mini,
  taxa_sums(data_fungi_mini) > 12000
))

test_that("is_mafft_installed returns a single logical", {
  res <- is_mafft_installed()
  expect_type(res, "logical")
  expect_length(res, 1)
  expect_false(is_mafft_installed(path = tempfile("no_mafft_here_")))
})

test_that("align_pq aligns a phyloseq refseq slot with DECIPHER", {
  skip_if_not_installed("DECIPHER")
  skip_if_not_installed("Biostrings")
  ali <- align_pq(df, method = "decipher")
  expect_s4_class(ali, "DNAStringSet")
  expect_length(ali, phyloseq::ntaxa(df))
  expect_length(unique(Biostrings::width(ali)), 1)
  expect_setequal(names(ali), phyloseq::taxa_names(df))
})

test_that("align_pq accepts a DNAStringSet and a DNAbin", {
  skip_if_not_installed("DECIPHER")
  dna <- Biostrings::DNAStringSet(phyloseq::refseq(df))
  expect_s4_class(align_pq(dna, method = "decipher"), "DNAStringSet")
  expect_s4_class(
    align_pq(ape::as.DNAbin(dna), method = "decipher"),
    "DNAStringSet"
  )
})

test_that("align_pq returns already-aligned sequences untouched", {
  aligned <- Biostrings::DNAStringSet(c(a = "AC-GT", b = "ACGGT"))
  expect_identical(align_pq(aligned), aligned)
  expect_message(align_pq(aligned, verbose = TRUE), "already share a width")
})

test_that("align_pq errors on bad input", {
  expect_error(align_pq(1:10), "phyloseq")
  expect_error(
    align_pq(Biostrings::DNAStringSet(c(a = "ACGT"))),
    "At least two sequences"
  )
  expect_error(
    align_pq(data_fungi_mini@otu_table),
    "phyloseq"
  )
})

test_that("align_pq aligns with mafft when it is installed", {
  skip_if_not(is_mafft_installed())
  skip_if_not_installed("ips")
  ali <- align_pq(df, method = "mafft")
  expect_s4_class(ali, "DNAStringSet")
  expect_length(ali, phyloseq::ntaxa(df))
  expect_length(unique(Biostrings::width(ali)), 1)
  expect_setequal(names(ali), phyloseq::taxa_names(df))
  expect_false(grepl("[a-z]", as.character(ali[[1]])))
})

test_that("align_pq reports a missing mafft executable", {
  skip_if_not_installed("ips")
  expect_error(
    align_pq(df, method = "mafft", exec = tempfile("no_mafft_here_")),
    "not found"
  )
})
