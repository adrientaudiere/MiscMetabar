test_that("fastq_to_fasta converts a plain FASTQ file", {
  fq <- system.file("extdata", "ex.fastq", package = "MiscMetabar")
  skip_if(fq == "")
  out_dir <- withr::local_tempdir()
  out <- fastq_to_fasta(fq, output_folder = out_dir)
  expect_true(file.exists(out))
  expect_match(out, "\\.fasta$")
  lines <- readLines(out)
  expect_true(grepl("^>", lines[1]))
  # one header + one sequence per FASTQ record
  n_records <- length(readLines(fq)) / 4
  expect_equal(length(lines), n_records * 2)
})

test_that("fastq_to_fasta converts a gzipped FASTQ file", {
  fq <- system.file("extdata", "ex_R1_001.fastq.gz", package = "MiscMetabar")
  skip_if(fq == "")
  out_dir <- withr::local_tempdir()
  out <- fastq_to_fasta(fq, output_folder = out_dir)
  expect_true(file.exists(out))
  expect_match(out, "ex_R1_001\\.fasta$")
  expect_true(grepl("^>", readLines(out, n = 1)))
})

test_that("fastq_to_fasta handles multiple files at once", {
  fq1 <- system.file("extdata", "ex.fastq", package = "MiscMetabar")
  fq2 <- system.file("extdata", "ex_R1_001.fastq.gz", package = "MiscMetabar")
  skip_if(fq1 == "" || fq2 == "")
  out_dir <- withr::local_tempdir()
  out <- fastq_to_fasta(c(fq1, fq2), output_folder = out_dir)
  expect_length(out, 2)
  expect_true(all(file.exists(out)))
})

test_that("fastq_to_fasta refuses to overwrite unless force = TRUE", {
  fq <- system.file("extdata", "ex.fastq", package = "MiscMetabar")
  skip_if(fq == "")
  out_dir <- withr::local_tempdir()
  fastq_to_fasta(fq, output_folder = out_dir)
  expect_error(fastq_to_fasta(fq, output_folder = out_dir), "already exists")
  expect_silent(suppressMessages(
    fastq_to_fasta(fq, output_folder = out_dir, force = TRUE)
  ))
})

test_that("fastq_to_fasta errors on missing input", {
  expect_error(
    fastq_to_fasta("does_not_exist.fastq", output_folder = tempdir()),
    "do not exist"
  )
})
