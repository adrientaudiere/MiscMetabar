data("enterotype")
sam_name_factice <- gsub("TS1_V2", "TS10_V2", sample_names(enterotype))


test_that("dist_pos_control function works fine", {
  expect_type(dist_pos_control(enterotype, sam_name_factice), "list")
  expect_s3_class(dist_pos_control(enterotype, sam_name_factice)[[1]], "data.frame")
  expect_s3_class(dist_pos_control(enterotype, sam_name_factice)[[2]], "data.frame")
  expect_equal(length(dist_pos_control(enterotype, sam_name_factice)), 2)
})

data(data_fungi)

res_seq <- suppressWarnings(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "cutoff_seq"
))
res_mixt <- suppressWarnings(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "cutoff_mixt"
))
res_diff <- suppressWarnings(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "cutoff_diff", min_diff_for_cutoff = 2
))
res_min <- suppressWarnings(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "min", min_diff_for_cutoff = 2
))
res_max <- suppressWarnings(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "max", min_diff_for_cutoff = 2
))
res_mean <- suppressWarnings(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "mean", min_diff_for_cutoff = 2
))


test_that("subset_taxa_tax_control function works fine", {
  expect_error(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]), method = "cutoff_diff"))
  expect_error(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]), method = "cut", min_diff_for_cutoff = 2))
  expect_s4_class(res_seq, "phyloseq")
  expect_s4_class(res_mixt, "phyloseq")
  expect_s4_class(res_diff, "phyloseq")
  expect_s4_class(res_min, "phyloseq")
  expect_s4_class(res_max, "phyloseq")
  expect_s4_class(res_mean, "phyloseq")
})

library(Biostrings)
test_that("search_exact_seq_pq function works fine", {
  expect_silent(search_primers <- search_exact_seq_pq(data_fungi, sequences = Biostrings::DNAStringSet(c("TTGAACGCACATTGCGCC", "ATCCCTACCTGATCCGAG"))))
  expect_equal(search_primers[[1]][3, 1], "932")
})
