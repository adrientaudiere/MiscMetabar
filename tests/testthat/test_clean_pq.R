data_fungi_test <- data_fungi
data_fungi_test@otu_table[, 1] <- rep(0, nrow(data_fungi_test@otu_table))
data_fungi_test@otu_table[10, ] <- rep(0, ncol(data_fungi_test@otu_table))

test_that("clean_pq return a phyloseq object after cleaning empty taxa and samples", {
  expect_s4_class(data_fungi_test, "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test), "phyloseq")
  skip_on_cran()
  expect_s4_class(clean_pq(data_fungi_test, verbose = TRUE), "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test, reorder_taxa = TRUE), "phyloseq")
  expect_s4_class(clean_pq(data_fungi_test, rename_taxa = TRUE), "phyloseq")
  expect_error(clean_pq(
    data_fungi_test,
    force_taxa_as_columns = TRUE,
    force_taxa_as_rows = TRUE
  ))
})

test_that("clean_pq clean empty taxa and samples", {
  expect_message(
    clean_pq(data_fungi_test),
    "Cleaning suppress 2 taxa and 1 samples."
  )
  skip_on_cran()
  expect_identical(
    nrow(data_fungi_test@otu_table) -
      nrow(clean_pq(data_fungi_test)@otu_table),
    1L
  )
  expect_identical(
    ncol(data_fungi_test@otu_table) -
      ncol(clean_pq(data_fungi_test)@otu_table),
    2L
  )
})

test_that("clean_pq force taxa in column", {
  expect_true(taxa_are_rows(clean_pq(data_fungi, force_taxa_as_rows = TRUE)))
  expect_false(taxa_are_rows(clean_pq(
    data_fungi,
    force_taxa_as_columns = TRUE
  )))
})

data_fungi_test2 <- data_fungi_test
taxa_names(data_fungi_test2) <- paste0(
  "ASV",
  seq_along(taxa_names(data_fungi_test2))
)
test_that("clean_pq works fine with bad taxa_names", {
  expect_s4_class(clean_pq(data_fungi_test2), "phyloseq")
  expect_message(
    clean_pq(data_fungi_test2),
    "Cleaning suppress 2 taxa and 1 samples"
  )
})

data_fungi_test3 <- data_fungi_test
taxa_names(data_fungi_test3@refseq) <- paste0(
  "ASV",
  seq_along(taxa_names(data_fungi_test3))
)
test_that("clean_pq works fine with bad taxa_names in refseq", {
  expect_s4_class(clean_pq(data_fungi_test3), "phyloseq")
  expect_message(
    clean_pq(data_fungi_test3),
    "Cleaning suppress 2 taxa and 1 samples"
  )
})

data_fungi_test4 <- data_fungi_test
taxa_names(data_fungi_test4@tax_table) <- paste0(
  "ASV",
  seq_along(taxa_names(data_fungi_test4))
)
test_that("clean_pq works fine with bad taxa_names in tax_table", {
  expect_s4_class(clean_pq(data_fungi_test4), "phyloseq")
  expect_message(
    clean_pq(data_fungi_test4),
    "Cleaning suppress 2 taxa and 1 samples"
  )
})

data_fungi_test5 <- data_fungi_test
sample_names(data_fungi_test5@tax_table) <- paste0(
  "SAMP",
  seq_along(taxa_names(data_fungi_test5))
)
test_that("clean_pq works fine with bad sample_names in sam_data", {
  expect_s4_class(clean_pq(data_fungi_test5), "phyloseq")
  expect_message(
    clean_pq(data_fungi_test5),
    "Cleaning suppress 2 taxa and 1 samples"
  )
})


data_fungi_test6 <- data_fungi_test
sample_names(data_fungi_test6)[1] <- paste0(
  "0",
  sample_names(data_fungi_test6)[1]
)
test_that("clean_pq send the good message with one sample with a 0 at the start of the name", {
  expect_s4_class(clean_pq(data_fungi_test6), "phyloseq")
  expect_message(expect_message(
    clean_pq(data_fungi_test6),
    "Cleaning suppress 2 taxa and 1 samples"
  ))
})

data_fungi_tax <- data_fungi_mini
data_fungi_tax@tax_table[1, "Genus"] <- " Russula "
data_fungi_tax@tax_table[1, "Species"] <- "Russula_sp"
data_fungi_tax@tax_table[2, "Species"] <- "uncultured"

test_that("clean_pq tax_table toggles are FALSE by default and leave tax_table untouched", {
  res <- suppressMessages(clean_pq(data_fungi_tax))
  expect_identical(
    as.character(res@tax_table[1, "Genus"]),
    " Russula "
  )
  expect_identical(
    as.character(res@tax_table[1, "Species"]),
    "Russula_sp"
  )
  expect_identical(
    as.character(res@tax_table[2, "Species"]),
    "uncultured"
  )
})

test_that("clean_pq tax_remove_border_spaces trims tax_table values", {
  res <- suppressMessages(
    clean_pq(data_fungi_tax, tax_remove_border_spaces = TRUE)
  )
  expect_identical(as.character(res@tax_table[1, "Genus"]), "Russula")
  expect_identical(as.character(res@tax_table[2, "Species"]), "uncultured")
})

test_that("clean_pq tax_replace_to_NA replaces NA-like patterns", {
  res <- suppressMessages(clean_pq(data_fungi_tax, tax_replace_to_NA = TRUE))
  expect_true(is.na(res@tax_table[2, "Species"]))
})

test_that("clean_pq tax_redundant_suffix drops redundant '_sp' tips", {
  pq <- data_fungi_mini
  pq@tax_table[1, "Genus"] <- "Russula"
  pq@tax_table[1, "Species"] <- "Russula_sp"
  res <- suppressMessages(clean_pq(pq, tax_redundant_suffix = TRUE))
  expect_true(is.na(res@tax_table[1, "Species"]))
})

test_that("clean_pq tax_redundant_suffix accepts a custom suffix", {
  pq <- data_fungi_mini
  pq@tax_table[1, "Genus"] <- "Russula"
  pq@tax_table[1, "Species"] <- "Russula_var"
  res <- suppressMessages(clean_pq(pq, tax_redundant_suffix = "_var"))
  expect_true(is.na(res@tax_table[1, "Species"]))
})
