skip_on_cran()

test_that("find_itsxpress() honours the MiscMetabar.itsxpresspath option", {
  withr::local_options(MiscMetabar.itsxpresspath = "/some/where/itsxpress")
  expect_equal(find_itsxpress(), "/some/where/itsxpress")
})

test_that("is_itsxpress_installed() returns a logical", {
  expect_type(is_itsxpress_installed(), "logical")
  expect_length(is_itsxpress_installed(), 1)
})

ext_path <- system.file("extdata", package = "MiscMetabar")

test_that("itsxpress_pq() builds paired-end unmerged commands for DADA2", {
  cmd <- itsxpress_pq(
    ext_path,
    folder_output = file.path(tempdir(), "itsxpress_test_pe"),
    region = "ITS2",
    taxa = "Fungi",
    cmd_is_run = FALSE
  )
  expect_type(cmd, "list")
  expect_length(cmd, 1)
  expect_match(cmd[[1]], "itsxpress")
  expect_match(cmd[[1]], "--fastq ")
  expect_match(cmd[[1]], "--fastq2 ")
  expect_match(cmd[[1]], "--outfile ")
  expect_match(cmd[[1]], "--outfile2 ")
  expect_match(cmd[[1]], "--region ITS2")
  expect_match(cmd[[1]], "--taxa Fungi")
  expect_false(grepl("--single_end", cmd[[1]]))
})

test_that("itsxpress_pq() supports merged paired-end and single-end modes", {
  cmd_merged <- itsxpress_pq(
    ext_path,
    folder_output = file.path(tempdir(), "itsxpress_test_merged"),
    single_file_output = TRUE,
    cmd_is_run = FALSE
  )
  expect_match(cmd_merged[[1]], "--fastq2 ")
  expect_false(grepl("--outfile2 ", cmd_merged[[1]]))

  cmd_se <- itsxpress_pq(
    ext_path,
    folder_output = file.path(tempdir(), "itsxpress_test_se"),
    paired_end = FALSE,
    pattern = "fastq",
    pattern_R1 = "ex.fastq",
    reversed_primers = TRUE,
    trim_ccs = TRUE,
    extra_itsxpress_args = "--keeptemp",
    cmd_is_run = FALSE
  )
  expect_match(cmd_se[[1]], "--single_end")
  expect_match(cmd_se[[1]], "--reversed_primers")
  expect_match(cmd_se[[1]], "--trim-ccs")
  expect_match(cmd_se[[1]], "--keeptemp")
  expect_false(grepl("--fastq2 ", cmd_se[[1]]))
})

test_that("itsxpress_pq() validates cluster_id and trim_ccs", {
  expect_error(
    itsxpress_pq(ext_path, cluster_id = 0.90, cmd_is_run = FALSE),
    "cluster_id"
  )
  expect_error(
    itsxpress_pq(
      ext_path,
      paired_end = TRUE,
      trim_ccs = TRUE,
      cmd_is_run = FALSE
    ),
    "trim_ccs"
  )
})
