data(enterotype)


test_that("write_pq function works fine with enterotype dataset", {
  testFolder <- tempdir()
  unlink(list.files(testFolder, full.names = TRUE), recursive = TRUE)
  expect_silent(write_pq(enterotype, path = testFolder, silent = TRUE))
  skip_on_cran()
  expect_silent(write_pq(
    enterotype,
    path = testFolder,
    silent = TRUE,
    sam_data_first = TRUE
  ))
  expect_silent(write_pq(
    enterotype,
    path = testFolder,
    silent = TRUE,
    write_sam_data = TRUE
  ))
  expect_silent(write_pq(enterotype, path = testFolder))
  expect_silent(write_pq(enterotype, one_file = TRUE, path = testFolder))
  expect_s4_class(read_pq(testFolder, taxa_are_rows = TRUE), "phyloseq")
  new_data_entero <- read_pq(testFolder, taxa_are_rows = TRUE)
  expect_identical(ntaxa(new_data_entero) - ntaxa(enterotype), 0L)
  expect_identical(nsamples(new_data_entero) - nsamples(enterotype), 0L)
})


test_that("write_pq function works fine with data_fungi dataset", {
  testFolder <- tempdir()
  unlink(list.files(testFolder, full.names = TRUE), recursive = TRUE)
  expect_silent(write_pq(data_fungi, path = testFolder, silent = TRUE))
  skip_on_cran()
  expect_silent(write_pq(
    data_fungi,
    path = testFolder,
    silent = TRUE,
    sam_data_first = TRUE
  ))
  expect_silent(write_pq(
    data_fungi,
    path = testFolder,
    silent = TRUE,
    write_sam_data = FALSE
  ))
  expect_silent(write_pq(
    data_fungi,
    path = testFolder,
    silent = TRUE,
    sam_data_first = TRUE,
    one_file = TRUE
  ))
  expect_silent(write_pq(
    data_fungi,
    path = testFolder,
    silent = TRUE,
    write_sam_data = FALSE,
    one_file = TRUE
  ))
  expect_silent(write_pq(data_fungi, path = testFolder))
  expect_silent(write_pq(data_fungi, one_file = TRUE, path = testFolder))
  expect_s4_class(read_pq(testFolder), "phyloseq")
  new_data_fungi <- read_pq(testFolder)
  expect_identical(ntaxa(new_data_fungi) - ntaxa(data_fungi), 0L)
  expect_identical(nsamples(new_data_fungi) - nsamples(data_fungi), 0L)
})

test_that("write_pq one_file keeps sample data rows equal to a taxon abundance", {
  testFolder <- file.path(tempdir(), "write_pq_one_sample")
  one_sample <- clean_pq(
    prune_samples(sample_names(data_fungi_mini)[1], data_fungi_mini),
    silent = TRUE
  )
  counts <- as.vector(one_sample@otu_table)
  one_sample@sam_data$same_as_count <- as.character(counts[counts > 0][1])

  for (sam_first in c(FALSE, TRUE)) {
    expect_silent(write_pq(
      one_sample,
      path = testFolder,
      one_file = TRUE,
      sam_data_first = sam_first,
      silent = TRUE
    ))
    all_in_one <- read.delim(
      file.path(testFolder, "ASV_table_allInOne.csv"),
      row.names = 1
    )
    expect_identical(
      nrow(all_in_one),
      ntaxa(one_sample) + ncol(one_sample@sam_data)
    )
  }
  unlink(testFolder, recursive = TRUE)
})

test_that("save_pq function works fine with data_fungi dataset", {
  skip_on_os("windows")
  skip_on_cran()
  testFolder <- tempdir()
  unlink(list.files(testFolder, full.names = TRUE), recursive = TRUE)
  expect_silent(save_pq(data_fungi, path = testFolder, silent = TRUE))
  expect_silent(save_pq(data_fungi, path = testFolder))
  expect_length(list.files(testFolder), 6)
})
