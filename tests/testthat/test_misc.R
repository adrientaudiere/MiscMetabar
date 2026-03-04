data(data_fungi)
data("enterotype")

test_that("dist_bycol works fine", {
  expect_length(
    dist_bycol(
      data_fungi@otu_table,
      as_binary_otu_table(data_fungi)@otu_table
    ),
    2
  )
  skip_on_cran()
  expect_error(length(dist_bycol(
    data_fungi@otu_table,
    enterotype@otu_table
  )))
})

test_that("all_object_size works fine", {
  expect_type(all_object_size(), "double")
})

test_that("diff_fct_diff_class works fine", {
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Sample_id,
      numeric_fonction = sum,
      na.rm = TRUE
    ),
    17852
  )
  skip_on_cran()
  expect_equal(
    round(
      diff_fct_diff_class(
        data_fungi@sam_data$Time,
        numeric_fonction = mean,
        na.rm = TRUE
      ),
      2
    ),
    5.80
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height == "Low",
      logical_method = "TRUE_if_one"
    ),
    TRUE
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height == "Low",
      logical_method = "NA_if_not_all_TRUE"
    ),
    NA
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height == "Low",
      logical_method = "FALSE_if_not_all_TRUE"
    ),
    FALSE
  )
  expect_equal(
    diff_fct_diff_class(
      data_fungi@sam_data$Height,
      character_method = "unique_or_na"
    ),
    NA
  )
  expect_equal(
    diff_fct_diff_class(
      c("IE", "IE"),
      character_method = "unique_or_na"
    ),
    "IE"
  )
  expect_equal(
    diff_fct_diff_class(
      c("IE", "IE", "TE", "TE"),
      character_method = "more_frequent"
    ),
    "IE"
  )
  expect_equal(
    diff_fct_diff_class(
      c("IE", "IE", "TE", "TE"),
      character_method = "more_frequent_without_equality"
    ),
    NA
  )
})


withr::local_envvar(
  R_USER_CACHE_DIR = tempfile(),
  .local_envir = teardown_env()
)

test_that("add_funguild_info works fine", {
  skip_on_cran()
  data_f <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 5000)
  expect_silent(
    data_f <- add_funguild_info(
      data_f,
      taxLevels = c(
        "Domain",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      )
    )
  )
  expect_equal(dim(data_f@tax_table)[2], 24)
})


test_that("are_modality_even_depth works fine", {
  expect_equal(
    are_modality_even_depth(data_fungi, "Time")$statistic[[1]],
    62.143
  )
  expect_equal(
    are_modality_even_depth(rarefy_even_depth(data_fungi), "Time")$p.value,
    1
  )
  expect_silent(are_modality_even_depth(data_fungi, "Height", boxplot = TRUE))
})


test_that("as_binary_otu_table works fine", {
  expect_s4_class(as_binary_otu_table(enterotype), "phyloseq")
  bin_pq <- as_binary_otu_table(data_fungi)
  expect_true(all(bin_pq@otu_table %in% c(0, 1)))
  skip_on_cran()
  bin_pq_5 <- as_binary_otu_table(data_fungi, min_number = 5)
  expect_true(all(bin_pq_5@otu_table %in% c(0, 1)))
  expect_true(sum(bin_pq_5@otu_table) <= sum(bin_pq@otu_table))
  expect_error(as_binary_otu_table("not_a_phyloseq"))
})


test_that("simplify_taxo works fine", {
  d_fm <- data_fungi_mini
  d_fm@tax_table[, "Species"] <- paste0(
    rep(
      c("s__", "s:"),
      ntaxa(d_fm) / 2
    ),
    d_fm@tax_table[, "Species"]
  )

  expect_s4_class(simplify_taxo(d_fm), "phyloseq")
  skip_on_cran()
  simplified <- simplify_taxo(d_fm)
  expect_false(any(grepl("s__", simplified@tax_table[, "Species"])))
  expect_false(any(grepl("s:", simplified@tax_table[, "Species"])))
  expect_s4_class(simplify_taxo(d_fm, remove_NA = TRUE), "phyloseq")
})


test_that("get_file_extension works fine", {
  expect_equal(get_file_extension("test.fasta"), "fasta")
  expect_equal(
    suppressWarnings(get_file_extension("test.fastq.gz")),
    c("fastq", "gz")
  )
  skip_on_cran()
  expect_warning(get_file_extension("test.file.fasta"))
  expect_error(get_file_extension("test_without_extension"))
})


test_that("perc works fine", {
  expect_equal(perc(0.5), 50)
  expect_equal(perc(1, 2), 50)
  skip_on_cran()
  expect_equal(perc(0.567, accuracy = 1), 56.7)
  expect_equal(perc(0.5, add_symbol = TRUE), "50%")
  expect_equal(perc(1, 4, add_symbol = TRUE), "25%")
})


test_that("funky_color works fine", {
  expect_type(funky_color(5), "character")
  expect_length(funky_color(10), 10)
})


test_that("fac2col works fine", {
  test_fac <- factor(c("A", "B", "A", "C"))
  result <- fac2col(test_fac)
  expect_type(result, "character")
  expect_length(result, 4)
  skip_on_cran()
  result_seed <- fac2col(test_fac, seed = 123)
  expect_type(result_seed, "character")
  test_fac_na <- factor(c("A", "B", NA, "C"))
  result_na <- fac2col(test_fac_na)
  expect_equal(result_na[3], "grey")
})


test_that("transp works fine", {
  expect_type(transp("red"), "character")
  expect_type(transp(c("red", "blue")), "character")
  expect_length(transp(c("red", "blue", "green")), 3)
  result <- transp("red", alpha = 0.5)
  expect_true(grepl("#", result))
})
