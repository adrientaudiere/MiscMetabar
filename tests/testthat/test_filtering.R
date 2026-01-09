data(data_fungi)
data(enterotype)

test_that("filt_taxa_pq works", {
  skip_on_cran()
  # Filter taxa based on a simple condition
  result <- filt_taxa_pq(data_fungi, taxa_sums(data_fungi) > 100)
  expect_s4_class(result, "phyloseq")
  expect_true(ntaxa(result) <= ntaxa(data_fungi))
})

test_that("filt_taxa_wo_NA works", {
  skip_on_cran()
  result <- filt_taxa_wo_NA(data_fungi, "Phylum")
  expect_s4_class(result, "phyloseq")
  # Check that no NAs remain in the Phylum column
  if (ntaxa(result) > 0) {
    expect_false(any(is.na(tax_table(result)[, "Phylum"])))
  }
})

test_that("distri_1_taxa works", {
  skip_on_cran()
  p <- distri_1_taxa(data_fungi, taxa_names(data_fungi)[1])
  expect_s3_class(p, "ggplot")
})
