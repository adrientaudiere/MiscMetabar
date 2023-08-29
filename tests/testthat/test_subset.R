data("data_fungi")
cond_taxa <- grepl("Endophyte", data_fungi@tax_table[, "Guild"])
names(cond_taxa) <- taxa_names(data_fungi)

test_that("subset_taxa_pq works fine", {
  expect_s4_class(subset_taxa_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota"), "phyloseq")
  expect_equal(ntaxa(subset_taxa_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota")), 1066)
  expect_s4_class(subset_taxa_pq(data_fungi, cond_taxa), "phyloseq")
  expect_equal(ntaxa(subset_taxa_pq(data_fungi, cond_taxa)), 128)
  expect_error(subset_taxa_pq(data_fungi, data_fungi@sam_data[["Height"]] == "Low"))
})

cond_samp <- grepl("A1", data_fungi@sam_data[["Sample_names"]])
test_that("subset_samples_pq works fine", {
  expect_s4_class(subset_samples_pq(data_fungi, data_fungi@sam_data[["Height"]] == "Low"), "phyloseq")
  expect_s4_class(subset_samples_pq(data_fungi, cond_samp), "phyloseq")
  expect_error(subset_samples_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota"))
})
