data(data_fungi)

test_that("format2dada2 works", {
  skip_on_cran()
  # Test with simple taxonomy
  tax <- c("k__Fungi", "p__Ascomycota", "c__Dothideomycetes")
  result <- format2dada2(tax)
  expect_type(result, "character")
})

test_that("format2dada2_species works", {
  skip_on_cran()
  # Test with species level taxonomy
  tax <- c("k__Fungi", "p__Ascomycota", "c__Dothideomycetes", "s__Species_name")
  result <- format2dada2_species(tax)
  expect_type(result, "character")
})

test_that("format2sintax works", {
  skip_on_cran()
  # Test with simple taxonomy
  tax <- c("Kingdom", "Phylum", "Class")
  result <- format2sintax(tax)
  expect_type(result, "character")
})
