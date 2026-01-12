data(data_fungi)

test_that("funguild_assign works", {
  skip_on_cran()
data_fungi_FUNGUILD <- funguild_assign(as.data.frame(tax_table(data_fungi)),
 db_funguild = db, tax_col = "Genus_species")
   expect_type(result, "list")
  expect_equal(ncol(data_fungi_FUNGUILD), 23)
})

test_that("formattable_pq works", {
  skip_on_cran()
  if (requireNamespace("formattable", quietly = TRUE)) {
   result <- formattable_pq(
    data_fungi,
    "Height",
    min_nb_seq_taxa = 10000,
    formattable_args = list("Phylum" = FALSE),
    log10trans = TRUE
  )
    expect_s3_class(result, "formattable")
  }
})
