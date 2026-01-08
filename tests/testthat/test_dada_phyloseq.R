data(data_fungi)
data(data_fungi_mini)
data(enterotype)


test_that("taxa_as_columns works fine", {
  result <- taxa_as_columns(data_fungi_mini)
  expect_s4_class(result, "phyloseq")
  expect_false(taxa_are_rows(result))
})


test_that("taxa_as_rows works fine", {
  result <- taxa_as_rows(data_fungi_mini)
  expect_s4_class(result, "phyloseq")
  expect_true(taxa_are_rows(result))
})


test_that("normalize_prop_pq works fine", {
  result <- normalize_prop_pq(data_fungi_mini)
  expect_s4_class(result, "phyloseq")
  skip_on_cran()
  result_no_log <- normalize_prop_pq(data_fungi_mini, base_log = NULL)
  expect_s4_class(result_no_log, "phyloseq")
  
  result_log10 <- normalize_prop_pq(data_fungi_mini, base_log = 10)
  expect_s4_class(result_log10, "phyloseq")
})


test_that("add_info_to_sam_data works fine", {
  result <- add_info_to_sam_data(data_fungi_mini)
  expect_s4_class(result, "phyloseq")
  expect_true("nb_seq" %in% sample_variables(result))
  expect_true("nb_otu" %in% sample_variables(result))
  
  skip_on_cran()
  # Test with custom dataframe
  new_df <- data.frame(
    variable_1 = runif(n = nsamples(data_fungi_mini), min = 1, max = 20)
  )
  rownames(new_df) <- sample_names(data_fungi_mini)
  result2 <- add_info_to_sam_data(data_fungi_mini, new_df)
  expect_true("variable_1" %in% sample_variables(result2))
  
  # Test with mismatched rownames
  new_df_bad <- data.frame(
    variable_1 = runif(n = nsamples(data_fungi_mini), min = 1, max = 20)
  )
  rownames(new_df_bad) <- paste0("BAD_", seq_len(nsamples(data_fungi_mini)))
  expect_error(add_info_to_sam_data(data_fungi_mini, new_df_bad))
})


test_that("taxa_only_in_one_level works fine", {
  data_fungi_mini_woNA4height <- subset_samples(
    data_fungi_mini,
    !is.na(data_fungi_mini@sam_data$Height)
  )
  expect_type(
    taxa_only_in_one_level(data_fungi_mini_woNA4height, "Height", "High"),
    "character"
  )
  skip_on_cran()
  expect_type(
    taxa_only_in_one_level(data_fungi_mini_woNA4height, "Height", "Low"),
    "character"
  )
  result <- taxa_only_in_one_level(data_fungi_mini_woNA4height, "Height", "High",
    min_nb_seq_taxa = 100
  )
  expect_type(result, "character")
})


test_that("physeq_or_string_to_dna works fine", {
  # Test with phyloseq
  result <- physeq_or_string_to_dna(data_fungi_mini)
  expect_s4_class(result, "DNAStringSet")
  
  # Test with character vector
  sequences_ex <- c(
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTAATAACGAATTCATTGAATCA",
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTT"
  )
  result2 <- physeq_or_string_to_dna(dna_seq = sequences_ex)
  expect_s4_class(result2, "DNAStringSet")
  expect_length(result2, 2)
  
  skip_on_cran()
  # Test errors
  expect_error(physeq_or_string_to_dna(physeq = data_fungi_mini, dna_seq = sequences_ex))
  expect_error(physeq_or_string_to_dna())
})


test_that("psmelt_samples_pq works fine", {
  result <- psmelt_samples_pq(data_fungi_mini)
  expect_s3_class(result, "tbl_df")
  expect_true("Hill_0" %in% colnames(result))
  expect_true("Hill_1" %in% colnames(result))
  expect_true("Hill_2" %in% colnames(result))
  
  skip_on_cran()
  # Test without hill numbers
  result2 <- psmelt_samples_pq(data_fungi_mini, hill_scales = NULL)
  expect_s3_class(result2, "tbl_df")
  expect_false("Hill_0" %in% colnames(result2))
  
  # Test with taxa_ranks
  result3 <- psmelt_samples_pq(data_fungi_mini, taxa_ranks = c("Class", "Family"))
  expect_s3_class(result3, "tbl_df")
})


test_that("rarefy_sample_count_by_modality works fine", {
  result <- rarefy_sample_count_by_modality(data_fungi_mini, "Height", rngseed = 123)
  expect_s4_class(result, "phyloseq")
  
  skip_on_cran()
  # Check that samples are evenly distributed
  result_table <- table(result@sam_data$Height)
  expect_equal(length(unique(result_table)), 1)
})


test_that("filt_taxa_wo_NA works fine", {
  result <- filt_taxa_wo_NA(data_fungi)
  expect_s4_class(result, "phyloseq")
  expect_true(ntaxa(result) <= ntaxa(data_fungi))
  
  skip_on_cran()
  # Test with specific ranks
  result2 <- filt_taxa_wo_NA(data_fungi, taxa_ranks = c(1:3))
  expect_s4_class(result2, "phyloseq")
  
  # Test with n_NA parameter
  result3 <- filt_taxa_wo_NA(data_fungi, n_NA = 1)
  expect_s4_class(result3, "phyloseq")
})


test_that("add_dna_to_phyloseq works fine", {
  # Create a physeq without refseq
  physeq_no_refseq <- phyloseq(
    otu_table(data_fungi_mini@otu_table),
    sample_data(data_fungi_mini@sam_data),
    tax_table(data_fungi_mini@tax_table)
  )
  
  skip_on_cran()
  # This will only work if taxa names are DNA sequences, which they might not be
  expect_error(add_dna_to_phyloseq(physeq_no_refseq))
})


test_that("reorder_taxa_pq works fine", {
  new_order <- rev(taxa_names(data_fungi_mini))
  result <- reorder_taxa_pq(data_fungi_mini, new_order)
  expect_s4_class(result, "phyloseq")
  expect_equal(taxa_names(result), new_order)
})
