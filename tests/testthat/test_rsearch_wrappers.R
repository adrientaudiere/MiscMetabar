skip_on_cran()
skip_if_not_installed("Rsearch")

test_that("wrappers abort with a clear message when Rsearch is missing", {
  with_mocked_bindings(
    requireNamespace = function(pkg, quietly) FALSE,
    .package = "base",
    {
      expect_error(plot_read_quality("x.fastq"), "Rsearch")
      expect_error(plot_ee_rate_dist("x.fastq"), "Rsearch")
      expect_error(vs_fastx_uniques("x.fasta"), "Rsearch")
      expect_error(vs_uchime_ref("x.fasta", "db.fasta"), "Rsearch")
    }
  )
})

test_that("plot_read_quality and plot_ee_rate_dist return ggplot objects", {
  fq <- system.file("extdata", "ex.fastq", package = "MiscMetabar")

  p1 <- plot_read_quality(fq, plot_title = FALSE)
  expect_true(inherits(p1, c("ggplot", "gtable", "ggExtraPlot")))

  p2 <- plot_ee_rate_dist(fq)
  expect_s3_class(p2, "ggplot")
})

test_that("vs_fastx_uniques dereplicates a fasta file", {
  skip_if_not(
    MiscMetabar:::is_vsearch_installed(),
    "vsearch not installed"
  )
  out <- tempfile(fileext = ".fasta")
  expect_invisible(
    vs_fastx_uniques(
      system.file("extdata", "ex_little.fasta", package = "MiscMetabar"),
      fastx_output = out
    )
  )
  expect_true(file.exists(out))
  dna <- Biostrings::readDNAStringSet(out)
  expect_gte(length(dna), 1)
  expect_true(all(grepl(";size=", names(dna))))
})

test_that("vs_uchime_ref splits chimeric and non-chimeric sequences", {
  skip_if_not(
    MiscMetabar:::is_vsearch_installed(),
    "vsearch not installed"
  )
  nc <- tempfile(fileext = ".fasta")
  ch <- tempfile(fileext = ".fasta")
  expect_invisible(
    vs_uchime_ref(
      system.file("extdata", "ex.fasta", package = "MiscMetabar"),
      database = system.file(
        "extdata",
        "100_sp_UNITE_sh_general_release_dynamic.fasta",
        package = "MiscMetabar"
      ),
      nonchimeras = nc,
      chimeras = ch
    )
  )
  expect_true(file.exists(nc))
  expect_true(file.exists(ch))
})
