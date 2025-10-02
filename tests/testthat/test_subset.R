data("data_fungi")
cond_taxa <- grepl("Endophyte", data_fungi@tax_table[, "Guild"])
names(cond_taxa) <- taxa_names(data_fungi)

test_that("subset_taxa_pq works fine", {
  expect_s4_class(subset_taxa_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota"), "phyloseq")
  skip_on_cran()
  expect_equal(ntaxa(subset_taxa_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota")), 1066)
  expect_s4_class(subset_taxa_pq(data_fungi, cond_taxa), "phyloseq")
  expect_equal(ntaxa(subset_taxa_pq(data_fungi, cond_taxa)), 128)
  expect_equal(ntaxa(subset_taxa_pq(data_fungi, cond_taxa, taxa_names_from_physeq = TRUE)), 128)
  expect_error(subset_taxa_pq(data_fungi, data_fungi@sam_data[["Height"]] == "Low"))
})

test_that("subset_taxa_random works fine", {
  # Test basic functionality
  expect_s4_class(subset_taxa_random(data_fungi, 20), "phyloseq")
  expect_equal(ntaxa(subset_taxa_random(data_fungi, 20, seed = 123)), 20)
  
  skip_on_cran()
  
  # Test reproducibility with seed
  set1 <- subset_taxa_random(data_fungi, 20, seed = 456)
  set2 <- subset_taxa_random(data_fungi, 20, seed = 456)
  expect_equal(taxa_names(set1), taxa_names(set2))
  
  # Test that different seeds produce different results (most of the time)
  set3 <- subset_taxa_random(data_fungi, 20, seed = 789)
  expect_false(all(taxa_names(set1) == taxa_names(set3)))
  
  # Test error handling
  expect_error(subset_taxa_random(data_fungi, 0))
  expect_error(subset_taxa_random(data_fungi, -5))
  expect_error(subset_taxa_random(data_fungi, ntaxa(data_fungi) + 1))
  expect_error(subset_taxa_random(data_fungi, 1.5))
  expect_error(subset_taxa_random("not_a_phyloseq", 20))
  
  # Test selecting all taxa
  all_taxa_subset <- subset_taxa_random(data_fungi, ntaxa(data_fungi), seed = 123)
  expect_equal(ntaxa(all_taxa_subset), ntaxa(data_fungi))
  
  # Test selecting just 1 taxon
  one_taxon <- subset_taxa_random(data_fungi, 1, seed = 123)
  expect_equal(ntaxa(one_taxon), 1)
})

cond_samp <- grepl("A1", data_fungi@sam_data[["Sample_names"]])
test_that("subset_samples_pq works fine", {
  expect_s4_class(subset_samples_pq(data_fungi, data_fungi@sam_data[["Height"]] == "Low"), "phyloseq")
  skip_on_cran()
  expect_s4_class(subset_samples_pq(data_fungi, cond_samp), "phyloseq")
  expect_error(subset_samples_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota"))
  data_fungi2 <- data_fungi
  data_fungi2@sam_data <- NULL
  expect_message(subset_samples_pq(data_fungi2, cond_samp))
  expect_s4_class(suppressMessages(subset_samples_pq(data_fungi2, cond_samp)), "phyloseq")
})
