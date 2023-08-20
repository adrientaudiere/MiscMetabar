test_that("list_fastq_files function works fine", {
   expect_s3_class(list_fastq_files("inst/extdata", paired_end =F, pattern_r1=""), "list")
   expect_equal(length(list_fastq_files("inst/extdata", paired_end =F, pattern_r1="")), 3)
   expect_equal(length(list_fastq_files("inst/extdata", paired_end =F, pattern_r1="", nb_files = 2)), 2)
   expect_s3_class(list_fastq_files("inst/extdata/"), "list")
   expect_equal(length(list_fastq_files("inst/extdata/")), 2)
})


data(data_fungi)
test_that("rename_samples_otu_table function works fine", {
    expect_s3_class(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T), "phyloseq")
    expect_equal(nrow(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T)), 185)
    expect_equal(ncol(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T)), 1420)
    expect_equal(sample_names(rename_samples_otu_table(data_fungi@otu_table, as.character(1:185), taxa_are_rows = T)), as.character(1:185))
})