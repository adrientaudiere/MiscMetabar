
data("data_fungi_sp_known")
test_that("Functions return object of class phyloseq ", {
  expect_s4_class(asv2otu(data_fungi_sp_known), "phyloseq")
  expect_s4_class(asv2otu(data_fungi_sp_known, method = "vsearch"), "phyloseq")
  expect_s4_class(lulu_phyloseq(data_fungi_sp_known)$new_physeq, "phyloseq")
})
