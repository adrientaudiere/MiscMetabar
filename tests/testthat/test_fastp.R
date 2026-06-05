test_that("is_fastp_installed returns a logical", {
  expect_type(is_fastp_installed(), "logical")
  expect_length(is_fastp_installed(), 1)
})

ext_path <- system.file("extdata", package = "MiscMetabar")

test_that("fastp builds paired-end commands without running them", {
  skip_on_cran()
  cmd <- fastp(
    ext_path,
    folder_output = file.path(tempdir(), "fastp_test"),
    cmd_is_run = FALSE
  )
  expect_type(cmd, "list")
  expect_length(cmd, 1)
  expect_match(cmd[[1]], "fastp")
  expect_match(cmd[[1]], "-i ")
  expect_match(cmd[[1]], "-I ")
  expect_match(cmd[[1]], "--detect_adapter_for_pe")
  expect_match(cmd[[1]], "--correction")
})

test_that("fastp single-end commands omit paired-end-only flags", {
  skip_on_cran()
  cmd <- fastp(
    ext_path,
    folder_output = file.path(tempdir(), "fastp_test_se"),
    paired_end = FALSE,
    pattern = "fastq",
    pattern_R1 = "ex.fastq",
    cmd_is_run = FALSE
  )
  expect_type(cmd, "list")
  expect_false(grepl("-I ", cmd[[1]]))
  expect_false(grepl("--detect_adapter_for_pe", cmd[[1]]))
})

test_that("fastp passes extra arguments through", {
  skip_on_cran()
  cmd <- fastp(
    ext_path,
    folder_output = file.path(tempdir(), "fastp_test_extra"),
    cmd_is_run = FALSE,
    extra_fastp_args = "--dedup"
  )
  expect_match(cmd[[1]], "--dedup")
})
