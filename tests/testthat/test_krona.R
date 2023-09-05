data("GlobalPatterns")
GA <- subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")

krona_error_or_not <- try(system("ktImportText 2>&1", intern = TRUE))

if (class(krona_error_or_not) == "try-error") {
  test_that("krona send an error when krona is not installed", {
    expect_warning(krona(GA, "Number.of.sequences.html"))
    expect_warning(krona(GA, "Number.of.ASVs.html", nb_seq = FALSE))
    expect_warning(merge_krona(c("Number.of.sequences.html", "Number.of.ASVs.html")))
  })
} else {
  test_that("krona function works fine with GlobalPatterns dataset", {
    testFolder <- tempdir()
    suppressWarnings(unlink(list.files(testFolder, full.names = TRUE), recursive = TRUE))
    expect_silent(krona(GA, file = paste0(testFolder, "/Number.of.sequences.html")))
    expect_silent(krona(GA, file = paste0(testFolder, "/Number.of.ASVs.html"), nb_seq = FALSE))
    expect_silent(merge_krona(
      c(
        paste0(testFolder, "/Number.of.sequences.html"),
        paste0(testFolder, "/Number.of.ASVs.html")
      ),
      output = paste0(testFolder, "mergeKrona.html")
    ))
    unlink(list.files(testFolder, full.names = TRUE), recursive = TRUE)
  })
}
