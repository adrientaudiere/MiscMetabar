data(data_fungi)
data(enterotype)

test_that("taxa_as_columns works", {
  result <- taxa_as_columns(data_fungi)
  expect_s4_class(result, "phyloseq")
  expect_false(taxa_are_rows(result))
})

test_that("taxa_as_rows works", {
  result <- taxa_as_rows(data_fungi)
  expect_s4_class(result, "phyloseq")
  expect_true(taxa_are_rows(result))
})

test_that("taxa_only_in_one_level works", {
  suppressWarnings(
    result <- taxa_only_in_one_level(data_fungi, "Height", "Low")
  )
  expect_length(result, 124)
})

test_that("tbl_sum_taxtable works", {
  if (requireNamespace("gtsummary", quietly = TRUE)) {
    result <- tbl_sum_taxtable(data_fungi)
    expect_s3_class(result, "tbl_summary")
  }
})

test_that("resolve_vector_ranks with method='preference' works and errors correctly", {
  expect_error(
    resolve_vector_ranks(c("a", "b", "a"), method = "preference"),
    "preference_index"
  )
  expect_error(
    resolve_vector_ranks(
      c("a", "b"),
      method = "preference",
      preference_index = 5
    ),
    "preference_index is higher"
  )
  expect_identical(
    resolve_vector_ranks(
      c("a", "b", "a"),
      method = "preference",
      preference_index = 1
    ),
    "a"
  )
  expect_identical(
    resolve_vector_ranks(
      c(NA, "b", "b"),
      method = "preference",
      preference_index = 1,
      second_method = "unanimity"
    ),
    "b"
  )
})

test_that("format2sintax works with taxnames", {
  names_in <- c(
    "seq1 k__Fungi;p__Basidiomycota;c__Agaricomycetes",
    "seq2 k__Fungi;p__Ascomycota;c__Sordariomycetes"
  )
  result <- format2sintax(taxnames = names_in)
  expect_type(result, "character")
  expect_length(result, 2)
  expect_true(all(grepl("tax=k:", result)))
})

test_that("format2dada2 errors when both/neither taxnames and fasta_db given", {
  expect_error(format2dada2(), "taxnames or fasta_db")
  expect_error(
    format2dada2(taxnames = "x", fasta_db = "y"),
    "not both"
  )
})

test_that("format2dada2 works with taxnames (from_sintax = TRUE)", {
  sintax_names <- c(
    "seq1;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales"
  )
  result <- format2dada2(taxnames = sintax_names, from_sintax = TRUE)
  expect_type(result, "character")
  expect_length(result, 1)
  expect_true(grepl("__", result))
})

test_that("format2dada2_species errors when both/neither given", {
  expect_error(format2dada2_species(), "taxnames or fasta_db")
  expect_error(
    format2dada2_species(taxnames = "x", fasta_db = "y"),
    "not both"
  )
})

test_that("resolve_vector_ranks replace_collapsed_rank_by_NA works", {
  expect_identical(
    resolve_vector_ranks(
      c("a", "b"),
      method = "consensus",
      replace_collapsed_rank_by_NA = TRUE
    ),
    NA_character_
  )
  expect_identical(
    resolve_vector_ranks(
      c("a", "a"),
      method = "consensus",
      replace_collapsed_rank_by_NA = TRUE
    ),
    "a"
  )
})
