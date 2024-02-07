suppressWarnings(cutadapt_error_or_not <-
  try(system("cutadapt 2>&1", intern = TRUE), silent = TRUE))

if (class(cutadapt_error_or_not) == "try-error" ||
    sum(grepl("ModuleNotFoundError", cutadapt_error_or_not)) != 0) {
  message(
    "cutadapt_remove_primers()  can't be
    tested when cutadapt is not installed"
  )
} else {
  test_that("cutadapt_remove_primers works fine", {
    expect_silent(suppressMessages(
      res_cut <- cutadapt_remove_primers("inst/extdata", "TTC", "GAA",
      folder_output = paste0(tempdir(), "/output_cutadapt")
    )))
    expect_type(res_cut, "list")
    expect_equal(length(res_cut), 1)
    expect_type(res_cut[[1]], "character")

    expect_equal(length(list.files(paste0(tempdir(), "/output_cutadapt"))), 2)

    expect_silent(
      suppressMessages(res_cut2 <- cutadapt_remove_primers(
        system.file("extdata",
                    package = "dada2"
        ),
        pattern_R1 = "F.fastq.gz",
        pattern_R2 = "R.fastq.gz",
        primer_fw = "TTC",
        primer_rev = "GAA",
        folder_output = tempdir()
      ))
    )

    expect_type(res_cut2, "list")
    expect_equal(length(res_cut2), 2)

    expect_type(res_cut2[[1]], "character")

    expect_equal(length(list.files(paste0(tempdir(), "/output_cutadapt"))), 2)

    expect_silent(
      suppressMessages(res_cut3 <- cutadapt_remove_primers(
        system.file("extdata",
                    package = "dada2"
        ),
        pattern_R1 = "F.fastq.gz",
        primer_fw = "TTC",
        folder_output = tempdir(),
        cmd_is_run = FALSE
      )
    ))

    expect_type(res_cut3, "list")
    expect_equal(length(res_cut3), 2)
    expect_type(res_cut3[[1]], "character")

    expect_equal(length(list.files(paste0(tempdir(), "/output_cutadapt"))), 2)
  })

  unlink(tempdir(), recursive = TRUE)

}
