data(enterotype)
data(data_fungi)

test_that("list_fastq_files function works fine", {
  expect_type(list_fastq_files("inst/extdata", paired_end = FALSE, pattern_R1 = ""), "list")
  expect_equal(length(unlist(list_fastq_files("inst/extdata", paired_end = FALSE, pattern_R1 = ""))), 3)
  expect_equal(length(unlist(list_fastq_files("inst/extdata", paired_end = FALSE, pattern_R1 = "", nb_files = 2))), 2)
  expect_type(list_fastq_files("inst/extdata/"), "list")
  expect_equal(length(list_fastq_files("inst/extdata/")), 2)
})

test_that("rename_samples_otu_table function works fine when taxa_are_rows", {
  expect_s4_class(rename_samples_otu_table(data_fungi, as.character(1:nsamples(data_fungi))), "otu_table")
  expect_equal(nrow(rename_samples_otu_table(data_fungi, as.character(1:nsamples(data_fungi)))), nsamples(data_fungi))
  expect_equal(ncol(rename_samples_otu_table(data_fungi, as.character(1:nsamples(data_fungi)))), ntaxa(data_fungi))
  expect_equal(sample_names(rename_samples_otu_table(data_fungi, as.character(1:nsamples(data_fungi)))), as.character(1:nsamples(data_fungi)))
  expect_error(rename_samples_otu_table(data_fungi, as.character(2:nsamples(data_fungi))))
})

data_fungi_row <- clean_pq(data_fungi, force_taxa_as_rows = TRUE)

test_that("rename_samples_otu_table function works fine when taxa_are_columns", {
  expect_s4_class(rename_samples_otu_table(data_fungi_row, as.character(1:nsamples(data_fungi_row))), "otu_table")
  expect_equal(ncol(rename_samples_otu_table(data_fungi_row, as.character(1:nsamples(data_fungi_row)))), nsamples(data_fungi_row))
  expect_equal(nrow(rename_samples_otu_table(data_fungi_row, as.character(1:nsamples(data_fungi_row)))), ntaxa(data_fungi_row))
  expect_equal(sample_names(rename_samples_otu_table(data_fungi_row, as.character(1:nsamples(data_fungi_row)))), as.character(1:nsamples(data_fungi)))
  expect_error(rename_samples_otu_table(data_fungi_row, as.character(2:nsamples(data_fungi_row))))
})

data_fungi_test <- data_fungi
data_fungi_test@otu_table[, 1] <- rep(0, nrow(data_fungi_test@otu_table))
data_fungi_test@otu_table[10, ] <- rep(0, ncol(data_fungi_test@otu_table))

test_that("track_wkflow function works fine", {
  expect_message(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype)))
  expect_s3_class(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype)), "data.frame")
  expect_s3_class(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype), obj_names = c("Fastq.files", "data_fungi", "Enterotype")), "data.frame")
  expect_s3_class(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, data_fungi_test), clean_pq = TRUE), "data.frame")
})

test_that("track_wkflow function works fine with taxonomy_rank", {
  expect_error(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype), taxonomy_rank = c(3, 5)))
  expect_s3_class(track_wkflow(list(data_fungi, enterotype), taxonomy_rank = c(3, 5)), "data.frame")
})

tree_A10_005 <- subset_samples(data_fungi, Tree_name == "A10-005")

test_that("track_wkflow_samples function works fine", {
  expect_message(track_wkflow_samples(tree_A10_005))
  expect_equal(length(track_wkflow_samples(tree_A10_005)), 3)
  expect_type(track_wkflow_samples(tree_A10_005), "list")
  expect_s3_class(track_wkflow_samples(tree_A10_005)[[1]], "data.frame")
})


derep_R1_001 <- dada2::derepFastq("inst/extdata/ex_R1_001.fastq.gz")
dada_R1_001 <- dada(derep_R1_001, selfConsist = TRUE, multithread = TRUE)
derep_R_001 <- dada2::derepFastq(c("inst/extdata/ex_R1_001.fastq.gz", "inst/extdata/ex_R2_001.fastq.gz"))
test_that("track_wkflow_samples function works fine with object of class matrix, dada and derep", {
  expect_s3_class(track_wkflow(list(data_fungi@otu_table, derep_R1_001, derep_R_001, "inst/extdata/ex_R1_001.fastq.gz", dada_R1_001)), "data.frame")
})

test_that("select_one_sample function works fine", {
  expect_message(A8_005 <- select_one_sample(data_fungi, "A8-005_S4_MERGED.fastq.gz"))
  expect_s4_class(A8_005, "phyloseq")
  expect_error(select_one_sample(data_fungi, "A8-005_S.fastq.gz"))
})
