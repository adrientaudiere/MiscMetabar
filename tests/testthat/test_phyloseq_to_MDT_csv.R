test_that("phyloseq_to_MDT_csv writes three CSV files in MDT orientation", {
  skip_on_cran()
  out_dir <- file.path(tempdir(), "mdt_csv_test")
  unlink(out_dir, recursive = TRUE)
  ps <- clean_pq(data_fungi_mini)
  files <- suppressMessages(
    phyloseq_to_MDT_csv(ps, path = out_dir, check_dwc = FALSE)
  )
  expect_named(files, c("OTU_table", "Taxonomy", "Samples"))
  expect_true(all(file.exists(files)))

  otu <- utils::read.csv(files[["OTU_table"]], check.names = FALSE)
  # OTU IDs in rows -> as many rows as taxa, leading id column
  expect_identical(colnames(otu)[1], "id")
  expect_identical(nrow(otu), ntaxa(ps))
  expect_identical(ncol(otu) - 1L, nsamples(ps))

  unlink(out_dir, recursive = TRUE)
})

test_that("phyloseq_to_MDT_csv writes TSV when sep = tab", {
  skip_on_cran()
  out_dir <- file.path(tempdir(), "mdt_tsv_test")
  unlink(out_dir, recursive = TRUE)
  files <- suppressMessages(
    phyloseq_to_MDT_csv(
      clean_pq(data_fungi_mini),
      path = out_dir,
      sep = "\t",
      check_dwc = FALSE
    )
  )
  expect_true(all(grepl("\\.tsv$", files)))
  expect_true(all(file.exists(files)))
  unlink(out_dir, recursive = TRUE)
})

test_that("phyloseq_to_MDT_csv warns about missing Darwin Core terms", {
  skip_on_cran()
  out_dir <- file.path(tempdir(), "mdt_dwc_test")
  unlink(out_dir, recursive = TRUE)
  expect_warning(
    phyloseq_to_MDT_csv(clean_pq(data_fungi_mini), path = out_dir),
    "Darwin Core"
  )
  unlink(out_dir, recursive = TRUE)
})
