# blast_to_phyloseq
# blast_pq
# filter_asv_blast
# add_blast_info

data("data_fungi")

blast_error_or_not <- try(system("blastn 2>&1", intern = TRUE))

if (class(blast_error_or_not) == "try-error") {
  test_that("Blast send an error when krona is not installed", {
    # expect_warning(XXX)
  })
} else {
  test_that("", 
  )
}