test_that("list_fastq_files function works fine", {
   expect_type(list_fastq_files("inst/extdata", paired_end =F, pattern_r1=""), "list")
   expect_equal(length(unlist(list_fastq_files("inst/extdata", paired_end =F, pattern_r1=""))), 3)
   expect_equal(length(unlist(list_fastq_files("inst/extdata", paired_end =F, pattern_r1="", nb_files = 2))), 2)
   expect_type(list_fastq_files("inst/extdata/"), "list")
   expect_equal(length(list_fastq_files("inst/extdata/")), 2)
})


data(data_fungi)
test_that("rename_samples_otu_table function works fine", {
    expect_s4_class(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T), "otu_table")
    expect_equal(nrow(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T)), 185)
    expect_equal(ncol(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T)), 1420)
    expect_equal(sample_names(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T)), as.character(1:185))
})

data(enterotype)
data_fungi_test <- data_fungi
data_fungi_test@otu_table[, 1] <- rep(0, nrow(data_fungi_test@otu_table))
data_fungi_test@otu_table[10, ] <- rep(0, ncol(data_fungi_test@otu_table))

test_that("track_wkflow function works fine", {
    expect_message(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype)))
    expect_s3_class(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype)), "data.frame")
    expect_s3_class(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype), obj_names = c("Fastq.files", "data_fungi", "Enterotype")), "data.frame")
    expect_s3_class(track_wkflow(list(unlist(list_fastq_files("inst/extdata/")), data_fungi, data_fungi_test), clean_pq = TRUE), "data.frame")
})

tree_A10_005 <- subset_samples(data_fungi, Tree_name == "A10-005")

test_that("track_wkflow function works fine", {
    expect_message(track_wkflow_samples(tree_A10_005))
    expect_equal(length(track_wkflow_samples(tree_A10_005)), 3)
    expect_type(track_wkflow_samples(tree_A10_005), "list")
    expect_s3_class(track_wkflow_samples(tree_A10_005)[[1]], "data.frame")
})
