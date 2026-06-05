test_that("phyloseq_to_MDT_excel writes a multi-sheet xlsx file", {
  skip_on_cran()
  skip_if_not_installed("writexl")
  out <- file.path(tempdir(), "mdt_test.xlsx")
  unlink(out)
  res <- suppressMessages(
    phyloseq_to_MDT_excel(
      clean_pq(data_fungi_mini),
      filename = out,
      check_dwc = FALSE
    )
  )
  expect_identical(res, out)
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
  unlink(out)
})

test_that("phyloseq_to_MDT_excel warns about missing Darwin Core terms", {
  skip_on_cran()
  skip_if_not_installed("writexl")
  out <- file.path(tempdir(), "mdt_dwc.xlsx")
  unlink(out)
  expect_warning(
    suppressMessages(
      phyloseq_to_MDT_excel(clean_pq(data_fungi_mini), filename = out)
    ),
    "decimalLatitude"
  )
  unlink(out)
})

test_that("phyloseq_to_MDT_excel check_dwc = FALSE is silent about DwC terms", {
  skip_on_cran()
  skip_if_not_installed("writexl")
  out <- file.path(tempdir(), "mdt_nodwc.xlsx")
  unlink(out)
  expect_no_warning(
    suppressMessages(
      phyloseq_to_MDT_excel(
        clean_pq(data_fungi_mini),
        filename = out,
        check_dwc = FALSE
      )
    )
  )
  unlink(out)
})
