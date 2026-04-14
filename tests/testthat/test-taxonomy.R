# .detect_tax_format --------------------------------------------------------

test_that(".detect_tax_format detects SINTAX format", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;tax=k:Fungi,p:Ascomycota,c:Sordariomycetes",
      "ATCGATCG",
      ">seq2;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes",
      "GCTAGCTA"
    ),
    tmp
  )
  expect_equal(MiscMetabar:::.detect_tax_format(tmp), "sintax")
  unlink(tmp)
})

test_that(".detect_tax_format detects UNITE format", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;k__Fungi;p__Ascomycota;c__Sordariomycetes",
      "ATCGATCG",
      ">seq2;k__Fungi;p__Basidiomycota;c__Agaricomycetes",
      "GCTAGCTA"
    ),
    tmp
  )
  expect_equal(MiscMetabar:::.detect_tax_format(tmp), "unite")
  unlink(tmp)
})

test_that(".detect_tax_format detects Greengenes2 format", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1 d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria",
      "ATCGATCG",
      ">seq2 d__Bacteria;p__Firmicutes;c__Bacilli",
      "GCTAGCTA"
    ),
    tmp
  )
  expect_equal(MiscMetabar:::.detect_tax_format(tmp), "greengenes2")
  unlink(tmp)
})

test_that(".detect_tax_format returns unknown for dada2 positional format", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">Fungi;Ascomycota;Sordariomycetes;Hypocreales;Nectriaceae;Fusarium;",
      "ATCGATCG",
      ">Fungi;Basidiomycota;Agaricomycetes;Agaricales;Amanitaceae;Amanita;",
      "GCTAGCTA"
    ),
    tmp
  )
  expect_equal(MiscMetabar:::.detect_tax_format(tmp), "unknown")
  unlink(tmp)
})

test_that(".detect_tax_format returns unknown for empty file", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(character(0), tmp)
  expect_equal(MiscMetabar:::.detect_tax_format(tmp), "unknown")
  unlink(tmp)
})

test_that(".detect_tax_format works with gzipped files", {
  tmp <- tempfile(fileext = ".fasta.gz")
  con <- gzfile(tmp, "w")
  writeLines(
    c(
      ">seq1;tax=k:Fungi,p:Ascomycota",
      "ATCGATCG"
    ),
    con
  )
  close(con)
  expect_equal(MiscMetabar:::.detect_tax_format(tmp), "sintax")
  unlink(tmp)
})

# .validate_ref_format ------------------------------------------------------

test_that(".validate_ref_format passes for correct SINTAX format", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;tax=k:Fungi,p:Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  expect_no_error(
    MiscMetabar:::.validate_ref_format(tmp, "sintax", "test_func")
  )
  unlink(tmp)
})

test_that(".validate_ref_format passes for unknown format (dada2 positional)", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">Fungi;Ascomycota;Sordariomycetes;",
      "ATCGATCG"
    ),
    tmp
  )
  expect_no_error(
    MiscMetabar:::.validate_ref_format(tmp, "dada2", "test_func")
  )
  expect_no_error(
    MiscMetabar:::.validate_ref_format(tmp, "sintax", "test_func")
  )
  unlink(tmp)
})

test_that(".validate_ref_format errors when SINTAX expected but UNITE detected", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;k__Fungi;p__Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    MiscMetabar:::.validate_ref_format(tmp, "sintax", "test_func"),
    "SINTAX"
  )
  unlink(tmp)
})

test_that(".validate_ref_format errors when dada2 expected but SINTAX detected", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;tax=k:Fungi,p:Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    MiscMetabar:::.validate_ref_format(tmp, "dada2", "test_func"),
    "dada2"
  )
  unlink(tmp)
})

test_that(".validate_ref_format errors when dada2 expected but Greengenes2 detected", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1 d__Bacteria;p__Proteobacteria",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    MiscMetabar:::.validate_ref_format(tmp, "dada2", "test_func"),
    "dada2"
  )
  unlink(tmp)
})

test_that(".validate_ref_format error message suggests dbpq conversion", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;k__Fungi;p__Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    MiscMetabar:::.validate_ref_format(tmp, "sintax", "test_func"),
    "dbpq::format2sintax"
  )
  expect_error(
    MiscMetabar:::.validate_ref_format(tmp, "dada2", "test_func"),
    "dbpq::format2dada2"
  )
  unlink(tmp)
})

# Format validation in assign_* functions -----------------------------------

test_that("assign_dada2 errors on SINTAX-formatted ref_fasta", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;tax=k:Fungi,p:Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    assign_dada2(
      seq2search = Biostrings::DNAStringSet(c(ASV1 = "ATCG")),
      ref_fasta = tmp
    ),
    "dada2"
  )
  unlink(tmp)
})

test_that("assign_sintax errors on UNITE-formatted ref_fasta", {
  skip_if_not(is_vsearch_installed())
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;k__Fungi;p__Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    assign_sintax(
      seq2search = Biostrings::DNAStringSet(c(ASV1 = "ATCG")),
      ref_fasta = tmp
    ),
    "SINTAX"
  )
  unlink(tmp)
})

test_that("assign_vsearch_lca errors on Greengenes2-formatted ref_fasta", {
  skip_if_not(is_vsearch_installed())
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1 d__Bacteria;p__Proteobacteria",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    assign_vsearch_lca(
      seq2search = Biostrings::DNAStringSet(c(ASV1 = "ATCG")),
      ref_fasta = tmp
    ),
    "SINTAX"
  )
  unlink(tmp)
})

test_that("assign_blastn errors on UNITE-formatted ref_fasta", {
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;k__Fungi;p__Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  data(data_fungi_mini)
  expect_error(
    assign_blastn(
      physeq = data_fungi_mini,
      ref_fasta = tmp
    ),
    "SINTAX"
  )
  unlink(tmp)
})

test_that("assign_mmseqs2 errors on UNITE-formatted ref_fasta", {
  skip_if_not(is_mmseqs2_installed())
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;k__Fungi;p__Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  expect_error(
    assign_mmseqs2(
      seq2search = Biostrings::DNAStringSet(c(ASV1 = "ATCG")),
      ref_fasta = tmp
    ),
    "SINTAX"
  )
  unlink(tmp)
})

test_that("assign_dada2 does not error when from_sintax = TRUE", {
  skip_on_cran()
  tmp <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">seq1;tax=k:Fungi,p:Ascomycota",
      "ATCGATCG"
    ),
    tmp
  )
  # Should not error on format validation (will fail later on actual assignment)
  # We just check the validation is skipped
  expect_error(
    assign_dada2(
      seq2search = Biostrings::DNAStringSet(c(ASV1 = "ATCG")),
      ref_fasta = tmp,
      from_sintax = TRUE
    ),
    regexp = "(?!dada2 format)", # not a format error
    perl = TRUE
  )
  unlink(tmp)
})
