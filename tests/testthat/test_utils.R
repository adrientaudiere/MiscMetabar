data(data_fungi)
data(enterotype)

test_that("unique_or_na works with default method", {
  expect_equal(unique_or_na(c("a", "a", "a")), "a")
  expect_true(is.na(unique_or_na(c("a", "b", "c"))))
  expect_equal(unique_or_na(c(1, 1, 1)), 1)
  expect_true(is.na(unique_or_na(c(1, 2, 3))))
})

test_that("unique_or_na works with factors", {
  f <- factor(c("a", "a"), ordered = TRUE)
  expect_equal(unique_or_na(f), f[1])
  f2 <- factor(c("a", "b", "c"), ordered = TRUE)
  result <- unique_or_na(f2)
  expect_true(is.na(result))
  expect_true(is.factor(result))
})

test_that("transp adds transparency to colors", {
  col <- transp("red", alpha = 0.5)
  expect_type(col, "character")
  expect_length(col, 1)
  
  cols <- transp(c("red", "blue"), alpha = 0.8)
  expect_length(cols, 2)
})

test_that("no_legend returns ggplot theme", {
  expect_s3_class(no_legend(), "theme")
})

test_that("funky_color generates color palette", {
  skip_on_cran()
  cols <- funky_color(5)
  expect_type(cols, "character")
  expect_length(cols, 5)
})

test_that("physeq_or_string_to_dna works", {
  skip_on_cran()
  result <- physeq_or_string_to_dna(data_fungi)
  expect_s4_class(result, "DNAStringSet")
})

test_that("resolve_vector_ranks works", {
  skip_on_cran()
  vec <- c("Kingdom", "Phylum", "Class")
  result <- resolve_vector_ranks(data_fungi, vec)
  expect_type(result, "integer")
})

test_that("sam_data_matching_names works", {
  skip_on_cran()
  result <- sam_data_matching_names(data_fungi)
  expect_true(is.logical(result))
})

test_that("rename_samples works", {
  skip_on_cran()
  new_names <- paste0("Sample_", seq_len(nsamples(data_fungi)))
  result <- rename_samples(data_fungi, new_names)
  expect_s4_class(result, "phyloseq")
  expect_equal(nsamples(result), nsamples(data_fungi))
})

test_that("is_falco_installed works", {
  result <- is_falco_installed()
  expect_type(result, "logical")
})

test_that("is_krona_installed works", {
  result <- is_krona_installed()
  expect_type(result, "logical")
})

test_that("is_mumu_installed works", {
  result <- is_mumu_installed()
  expect_type(result, "logical")
})
