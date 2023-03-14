data("data_fungi_sp_known")
test_that("Functions return object of class phyloseq ", {
  expect_s4_class(data_fungi_sp_known, "phyloseq")
  expect_s4_class(asv2otu(data_fungi_sp_known), "phyloseq")
  expect_s4_class(asv2otu(data_fungi_sp_known, method = "vsearch"), "phyloseq")
  expect_s4_class(lulu_pq(data_fungi_sp_known)$new_physeq, "phyloseq")
  expect_s4_class(as_binary_otu_table(data_fungi_sp_known), "phyloseq")
})
