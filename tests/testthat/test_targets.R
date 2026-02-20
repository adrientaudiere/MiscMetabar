data(enterotype)
data(data_fungi)

test_that("list_fastq_files function works fine", {
  skip_on_cran()
  expect_type(
    list_fastq_files(
      "inst/extdata",
      paired_end = FALSE,
      pattern_R1 = ""
    ),
    "list"
  )
  expect_equal(
    length(unlist(
      list_fastq_files(
        "inst/extdata",
        paired_end = FALSE,
        pattern_R1 = ""
      )
    )),
    3
  )
  expect_equal(
    length(unlist(
      list_fastq_files(
        "inst/extdata",
        paired_end = FALSE,
        pattern_R1 = "",
        nb_files = 2
      )
    )),
    2
  )
  expect_type(list_fastq_files("inst/extdata/"), "list")
  expect_equal(length(list_fastq_files("inst/extdata/")), 2)
})

test_that("rename_samples_otu_table function works fine when taxa_are_rows", {
  expect_s4_class(
    rename_samples_otu_table(data_fungi, as.character(1:nsamples(data_fungi))),
    "otu_table"
  )
  expect_equal(
    nrow(rename_samples_otu_table(
      data_fungi,
      as.character(
        1:nsamples(data_fungi)
      )
    )),
    nsamples(data_fungi)
  )
  expect_equal(
    ncol(rename_samples_otu_table(
      data_fungi,
      as.character(
        1:nsamples(data_fungi)
      )
    )),
    ntaxa(data_fungi)
  )
  expect_equal(
    sample_names(rename_samples_otu_table(
      data_fungi,
      as.character(
        1:nsamples(data_fungi)
      )
    )),
    as.character(1:nsamples(data_fungi))
  )
  expect_error(rename_samples_otu_table(
    data_fungi,
    as.character(2:nsamples(data_fungi))
  ))
})

data_fungi_row <- clean_pq(data_fungi, force_taxa_as_rows = TRUE)

test_that("rename_samples_otu_table function works fine when taxa_are_columns", {
  skip_on_cran()
  expect_s4_class(
    rename_samples_otu_table(
      data_fungi_row,
      as.character(1:nsamples(data_fungi_row))
    ),
    "otu_table"
  )
  expect_equal(
    ncol(rename_samples_otu_table(
      data_fungi_row,
      as.character(1:nsamples(data_fungi_row))
    )),
    nsamples(data_fungi_row)
  )
  expect_equal(
    nrow(rename_samples_otu_table(
      data_fungi_row,
      as.character(1:nsamples(data_fungi_row))
    )),
    ntaxa(data_fungi_row)
  )
  expect_equal(
    sample_names(rename_samples_otu_table(
      data_fungi_row,
      as.character(1:nsamples(data_fungi_row))
    )),
    as.character(1:nsamples(data_fungi))
  )
  expect_error(rename_samples_otu_table(
    data_fungi_row,
    as.character(2:nsamples(data_fungi_row))
  ))
})

data_fungi_test <- data_fungi
data_fungi_test@otu_table[, 1] <-
  rep(0, nrow(data_fungi_test@otu_table))
data_fungi_test@otu_table[10, ] <-
  rep(0, ncol(data_fungi_test@otu_table))

test_that("track_wkflow function works fine", {
  skip_on_os("windows")
  skip_on_cran()
  expect_message(track_wkflow(list(
    unlist(list_fastq_files("inst/extdata/")),
    data_fungi,
    enterotype
  )))
  expect_s3_class(
    track_wkflow(list(
      unlist(list_fastq_files("inst/extdata/")),
      data_fungi,
      enterotype
    )),
    "data.frame"
  )
  expect_s3_class(
    track_wkflow(
      list(unlist(list_fastq_files("inst/extdata/")), data_fungi, enterotype),
      obj_names = c("Fastq.files", "data_fungi", "Enterotype")
    ),
    "data.frame"
  )
  expect_s3_class(
    track_wkflow(
      list(
        unlist(list_fastq_files("inst/extdata/")),
        data_fungi,
        data_fungi_test
      ),
      clean_pq = TRUE
    ),
    "data.frame"
  )
})

test_that("track_wkflow function works fine with taxonomy_rank", {
  skip_on_cran()
  expect_error(track_wkflow(
    list(
      unlist(list_fastq_files("inst/extdata/")),
      data_fungi,
      enterotype
    ),
    taxonomy_rank = c(3, 5)
  ))
  expect_s3_class(
    track_wkflow(list(data_fungi, enterotype), taxonomy_rank = c(3, 5)),
    "data.frame"
  )
})

tree_A10_005 <- subset_samples(data_fungi, Tree_name == "A10-005")

test_that("track_wkflow_samples function works fine", {
  skip_on_cran()
  expect_message(track_wkflow_samples(tree_A10_005))
  expect_equal(length(track_wkflow_samples(tree_A10_005)), 3)
  expect_type(track_wkflow_samples(tree_A10_005), "list")
  expect_s3_class(track_wkflow_samples(tree_A10_005)[[1]], "data.frame")
})


derep_R1_001 <- dada2::derepFastq("inst/extdata/ex_R1_001.fastq.gz")
dada_R1_001 <-
  dada(derep_R1_001, selfConsist = TRUE, multithread = TRUE)
derep_R_001 <-
  dada2::derepFastq(c(
    "inst/extdata/ex_R1_001.fastq.gz",
    "inst/extdata/ex_R2_001.fastq.gz"
  ))
test_that("track_wkflow_samples function works fine with object of class matrix, dada and derep", {
  skip_on_os("windows")
  if (requireNamespace("pbapply")) {
    expect_s3_class(
      track_wkflow(
        list(
          data_fungi@otu_table,
          derep_R1_001,
          derep_R_001,
          "inst/extdata/ex_R1_001.fastq.gz",
          dada_R1_001
        )
      ),
      "data.frame"
    )
  }
})

test_that("track_wkflow_samples works with matrix input", {
  skip_on_cran()
  skip_on_os("windows")
  mat <- as(data_fungi@otu_table, "matrix")
  if (!taxa_are_rows(data_fungi)) {
    # rows are samples already
  } else {
    mat <- t(mat)
  }
  res <- track_wkflow_samples(mat)
  expect_type(res, "list")
  expect_equal(length(res), nrow(mat))
  expect_s3_class(res[[1]], "data.frame")
})

test_that("track_wkflow_samples works with single dada-class", {
  skip_on_cran()
  skip_on_os("windows")
  res <- track_wkflow_samples(list("sample1" = dada_R1_001))
  expect_type(res, "list")
  expect_equal(length(res), 1)
  expect_equal(names(res), "sample1")
  expect_s3_class(res[[1]], "data.frame")
})

test_that("track_wkflow_samples works with list of dada-class", {
  skip_on_cran()
  skip_on_os("windows")
  dada_list <- list("sampleA" = dada_R1_001, "sampleB" = dada_R1_001)
  res <- track_wkflow_samples(dada_list)
  expect_type(res, "list")
  expect_equal(length(res), 2)
  expect_true(all(c("sampleA", "sampleB") %in% names(res)))
  expect_s3_class(res[[1]], "data.frame")
})

test_that("track_wkflow_samples works with single derep-class", {
  skip_on_cran()
  skip_on_os("windows")
  res <- track_wkflow_samples(list("sample1" = derep_R1_001))
  expect_type(res, "list")
  expect_equal(length(res), 1)
  expect_equal(names(res), "sample1")
  expect_s3_class(res[[1]], "data.frame")
})

test_that("track_wkflow_samples works with list of derep-class", {
  skip_on_cran()
  skip_on_os("windows")
  res <- track_wkflow_samples(derep_R_001)
  expect_type(res, "list")
  expect_equal(length(res), length(derep_R_001))
  expect_s3_class(res[[1]], "data.frame")
})

test_that("track_wkflow_samples works with character vector (fastq paths)", {
  skip_on_cran()
  skip_on_os("windows")
  fastq_paths <- c(
    "inst/extdata/ex_R1_001.fastq.gz",
    "inst/extdata/ex_R2_001.fastq.gz"
  )
  res <- track_wkflow_samples(fastq_paths)
  expect_type(res, "list")
  expect_equal(length(res), 2)
  expect_equal(names(res), basename(fastq_paths))
  expect_s3_class(res[[1]], "data.frame")
})

test_that("track_wkflow_samples works with named character vector", {
  skip_on_cran()
  skip_on_os("windows")
  fastq_paths <- c(
    "sampleX" = "inst/extdata/ex_R1_001.fastq.gz",
    "sampleY" = "inst/extdata/ex_R2_001.fastq.gz"
  )
  res <- track_wkflow_samples(fastq_paths)
  expect_type(res, "list")
  expect_equal(length(res), 2)
  expect_equal(names(res), c("sampleX", "sampleY"))
})

test_that("track_wkflow_samples works with mixed object types", {
  skip_on_cran()
  skip_on_os("windows")
  mat <- as(data_fungi@otu_table, "matrix")
  if (taxa_are_rows(data_fungi)) {
    mat <- t(mat)
  }
  res <- track_wkflow_samples(list(
    tree_A10_005,
    mat
  ))
  expect_type(res, "list")
  expect_true(length(res) > 0)
  expect_s3_class(res[[1]], "data.frame")
})

test_that("select_one_sample function works fine", {
  expect_message(
    A8_005 <-
      select_one_sample(data_fungi, "A8-005_S4_MERGED.fastq.gz")
  )
  expect_s4_class(A8_005, "phyloseq")
  expect_error(select_one_sample(data_fungi, "A8-005_S.fastq.gz"))
})

test_that("subsample_fastq function works fine", {
  expect_silent(subsample_fastq(
    "inst/extdata/ex_R1_001.fastq.gz",
    "your_path_to_output"
  ))
  file.exists("your_path_to_output/ex_R1_001.fastq.gz")
  unlink("your_path_to_output", recursive = TRUE)
  expect_silent(subsample_fastq(
    list_fastq_files("inst/extdata"),
    "your_path_to_output2",
    nb_seq = 10
  ))
  file.exists("your_path_to_output2/ex_R1_001.fastq.gz")
  unlink("your_path_to_output2", recursive = TRUE)
})

test_that("sample_data_with_new_names function works fine", {
  sam_file <- system.file("extdata", "sam_data.csv", package = "MiscMetabar")
  expect_silent(
    newdf <- sample_data_with_new_names(
      sam_file,
      paste0("Samples_", seq(1, 185))
    )
  )
  expect_equal(dim(newdf)[1], 185)
  expect_equal(dim(newdf)[2], 7)
})


test_that("sample_data_with_new_names function works fine", {
  testFastqs_fw <- c(
    system.file("extdata", "sam1F.fastq.gz", package = "dada2"),
    system.file("extdata", "sam2F.fastq.gz", package = "dada2")
  )
  testFastqs_rev <- c(
    system.file("extdata", "sam1R.fastq.gz", package = "dada2"),
    system.file("extdata", "sam2R.fastq.gz", package = "dada2")
  )
  expect_silent(
    filt_fastq_fw <- filter_trim(testFastqs_fw, output_fw = tempdir())
  )
  expect_equal(length(derepFastq(filt_fastq_fw[1])), 2)
  expect_message(
    filt_fastq_pe <- filter_trim(
      fw = testFastqs_fw,
      rev = testFastqs_rev,
      output_fw = paste0(tempdir(), "/", "fw"),
      output_rev = paste0(tempdir(), "/", "rev")
    )
  )
  unlink(paste0(tempdir(), "/", "rev"))
  unlink(paste0(tempdir(), "/", "fw"))
  expect_equal(length(derepFastq(filt_fastq_pe[[1]])), 2)
  expect_equal(length(derepFastq(filt_fastq_pe[[2]])), 2)
})

test_that("add_info_to_sam_data function works fine with data_fungi", {
  new_df <- data.frame(
    variable_1 = runif(n = nsamples(data_fungi), min = 1, max = 20),
    variable_2 = runif(n = nsamples(data_fungi), min = 1, max = 2)
  )
  rownames(new_df) <- sample_names(data_fungi)
  expect_silent(data_fungi2 <- add_info_to_sam_data(data_fungi, new_df))
  expect_equal(dim(data_fungi2@sam_data)[2], 11)
  expect_equal(length(data_fungi2@sam_data$nb_seq), 185)
  expect_equal(length(data_fungi2@sam_data$nb_otu), 185)
})

unlink(paste0(tempdir(), "/", "rev"))
unlink(paste0(tempdir(), "/", "fw"))

test_that("rename_samples works with phyloseq object", {
  new_names <- paste0("S", seq_len(nsamples(data_fungi)))
  res <- rename_samples(data_fungi, new_names)
  expect_s4_class(res, "phyloseq")
  expect_equal(sample_names(res), new_names)
  expect_error(rename_samples(data_fungi, new_names[-1]))
})

test_that("rename_samples works with otu_table", {
  otu <- data_fungi@otu_table
  new_names <- paste0("S", seq_len(nsamples(data_fungi)))
  res <- rename_samples(otu, new_names)
  expect_s4_class(res, "otu_table")
  expect_equal(sample_names(res), new_names)
  expect_error(rename_samples(otu, new_names[-1]))
})

test_that("sam_data_matching_names returns expected structure", {
  sam_file <- system.file("extdata", "sam_data.csv", package = "MiscMetabar")
  extdata_path <- system.file("extdata", package = "MiscMetabar")
  res <- suppressWarnings(suppressMessages(sam_data_matching_names(
    sam_file,
    sample_col_name = "X",
    path_raw_seq = extdata_path,
    verbose = FALSE,
    sep = "\t"
  )))
  expect_type(res, "list")
  expect_true("sam_names_matching" %in% names(res))
  expect_true("sam_data" %in% names(res))
  expect_s3_class(res$sam_names_matching, "tbl_df")
})
