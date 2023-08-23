data("enterotype")
sam_name_factice <- gsub("TS1_V2", "TS10_V2", sample_names(enterotype))


test_that("dist_pos_control function works fine", {
  expect_type(dist_pos_control(enterotype, sam_name_factice), "list")
  expect_equal(length(dist_pos_control(enterotype, sam_name_factice)), 2)
})

data(data_fungi)

res_seq <- subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "cutoff_seq"
)
res_mixt <- subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "cutoff_mixt"
)
res_diff <- subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "cutoff_diff", min_diff_for_cutoff = 2
)
res_min <- subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "min", min_diff_for_cutoff = 2
)
res_max <- subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "max", min_diff_for_cutoff = 2
)
res_mean <- subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]),
  method = "mean", min_diff_for_cutoff = 2
)


test_that("subset_taxa_tax_control function works fine", {
  expect_error(subset_taxa_tax_control(data_fungi, as.numeric(data_fungi@otu_table[, 50]), method = "cutoff_diff"))
  expect_s4_class(res_seq, "phyloseq")
  expect_s4_class(res_mixt, "phyloseq")
  expect_s4_class(res_diff, "phyloseq")
  expect_s4_class(res_min, "phyloseq")
  expect_s4_class(res_max, "phyloseq")
  expect_s4_class(res_mean, "phyloseq")
})
