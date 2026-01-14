data(data_fungi)

test_that("format2dada2 works works with Unite", {
 result <- format2dada2(
    test_path(
      "inst/extdata",
      "sh_general_release_dynamic_19.02.2025_MINI.fasta"
    ),
    output_path = tempfile()
  )
  expect_s4_class(result, "DNAStringSet")

  expect_equal(names(result)[[1]], "Abrothallus_subhalei|MT153946|SH1227328.10FU|refs|k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Abrothallales;f__Abrothallaceae;g__Abrothallus;s__Abrothallus_subhalei;")
})


test_that("format2dada2 works works with Eukaryome", {
 result <- format2dada2(
    test_path(
      "inst/extdata",
      "General_EUK_ITS_v2.0_MINI.fasta"
    ),
    output_path = tempfile()
  )
  expect_s4_class(result, "DNAStringSet")

  expect_equal(names(result)[[1]], "AABX02000063;k__Straminipila;p__Oomycota;c__Peronosporomycetes;o__Peronosporales;f__Peronosporaceae;g__Phytophthora;s__unclassified;")
})


test_that("format2dada2_species works with Unite", {
   result <- format2dada2_species(
    test_path(
      "inst/extdata",
      "sh_general_release_dynamic_19.02.2025_MINI.fasta"
    ),
    output_path = tempfile()
  )
  expect_s4_class(result, "DNAStringSet")

  expect_equal(names(result)[[1]], "Abrothallus_subhalei|MT153946|SH1227328.10FU|refs| Abrothallus Abrothallus_subhalei")
})

test_that("format2dada2_species works with Eukaryome", {
   result <- format2dada2_species(
    test_path(
      "inst/extdata",
      "General_EUK_ITS_v2.0_MINI.fasta"
    ),
    output_path = tempfile()
  )
  expect_s4_class(result, "DNAStringSet")

  expect_equal(names(result)[[1]], "AABX02000063; Phytophthora unclassified")
})


test_that("format2sintax works with Unite", {
 result <- format2sintax(
    test_path(
      "inst/extdata",
      "sh_general_release_dynamic_19.02.2025_MINI.fasta"
    ),
    output_path = tempfile()
  )

  expect_s4_class(result, "DNAStringSet")
  expect_equal(names(result)[[1]], "Abrothallus_subhalei|MT153946|SH1227328.10FU|refs|;tax=k:Fungi,p:Ascomycota,c:Dothideomycetes,o:Abrothallales,f:Abrothallaceae,g:Abrothallus,s:Abrothallus_subhalei")

})

test_that("format2sintax works with Eukaryome", {
  result <- format2sintax(
    test_path(
      "inst/extdata",
      "General_EUK_ITS_v2.0_MINI.fasta"
    ),
    output_path = tempfile()
  )

  expect_s4_class(result, "DNAStringSet")
  expect_equal(names(result)[[1]], "AABX02000063,;tax=k:Straminipila,p:Oomycota,c:Peronosporomycetes,o:Peronosporales,f:Peronosporaceae,g:Phytophthora,s:unclassified")
})
