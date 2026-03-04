data(data_fungi)

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
  expect_true(ggplot2::is_theme(no_legend()[[1]]))
})

test_that("funky_color generates color palette", {
  cols <- funky_color(5)
  expect_type(cols, "character")
  expect_length(cols, 5)
  expect_error(funky_color(NA))
  expect_length(funky_color(0), 0)
})

test_that("physeq_or_string_to_dna works", {
  result <- physeq_or_string_to_dna(data_fungi)
  expect_s4_class(result, "DNAStringSet")
  sequences_ex <- c(
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTAATAACGAATTCATTGAATCA",
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTT",
    "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAGAGGTG",
    "TACCTATGTTGCCTTGGCGGCTAAACCTACC",
    "CGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCATGAGGTG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
    "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG"
  )
  dna2 <- physeq_or_string_to_dna(dna_seq = sequences_ex)

  expect_s4_class(dna2, "DNAStringSet")
})

test_that("resolve_vector_ranks works with a unique value", {
  expect_equal(resolve_vector_ranks(c("A")), "A")
  expect_equal(
    resolve_vector_ranks(c("A"), method = "preference", preference_index = 1),
    "A"
  )
  expect_equal(resolve_vector_ranks(c("A"), method = "abs_majority"), "A")
  expect_equal(resolve_vector_ranks(c("A"), method = "rel_majority"), "A")

  expect_equal(
    resolve_vector_ranks("A", method = "abs_majority", nb_agree_threshold = 0),
    "A"
  )
  expect_true(is.na(resolve_vector_ranks(
    "A",
    method = "abs_majority",
    nb_agree_threshold = 2
  )))

  expect_equal(
    resolve_vector_ranks("A", method = "rel_majority", nb_agree_threshold = 0),
    "A"
  )
  expect_equal(
    resolve_vector_ranks("A", method = "rel_majority", nb_agree_threshold = 1),
    "A"
  )
  expect_true(is.na(resolve_vector_ranks(
    "A",
    method = "rel_majority",
    nb_agree_threshold = 2
  )))

  expect_equal(resolve_vector_ranks(c("A"), method = "unanimity"), "A")
})


test_that("resolve_vector_ranks works with a vector of unique value", {
  vec <- c("A", "A", "A")
  expect_equal(resolve_vector_ranks(vec), "A")
  expect_equal(
    resolve_vector_ranks(
      vec,
      method = "preference",
      preference_index = 1
    ),
    "A"
  )
  expect_equal(
    resolve_vector_ranks(vec, method = "abs_majority"),
    "A"
  )
  expect_equal(
    resolve_vector_ranks(vec, method = "rel_majority"),
    "A"
  )
  expect_equal(
    resolve_vector_ranks(vec, method = "unanimity"),
    "A"
  )
})


test_that("resolve_vector_ranks works with a vector of NA", {
  vec <- c(NA, NA, NA)
  expect_true(is.na(resolve_vector_ranks(vec)))
  expect_true(is.na(resolve_vector_ranks(
    vec,
    method = "preference",
    preference_index = 1
  )))
  expect_true(is.na(resolve_vector_ranks(vec, method = "abs_majority")))
  expect_true(is.na(resolve_vector_ranks(vec, method = "rel_majority")))
  expect_true(is.na(resolve_vector_ranks(vec, method = "unanimity")))
})


test_that("resolve_vector_ranks works with a vector of 2 A and one NA", {
  vec <- c("A", "A", NA)
  expect_equal(resolve_vector_ranks(vec), "A")
  expect_equal(
    resolve_vector_ranks(vec, method = "preference", preference_index = 1),
    "A"
  )
  expect_equal(resolve_vector_ranks(vec, method = "abs_majority"), "A")
  expect_equal(resolve_vector_ranks(vec, method = "rel_majority"), "A")
  expect_equal(resolve_vector_ranks(vec, method = "unanimity"), "A")
  expect_true(is.na(resolve_vector_ranks(
    vec,
    method = "unanimity",
    strict = TRUE
  )))
})


test_that("resolve_vector_ranks works with a vector of one A, one B and one NA", {
  vec <- c("A", "B", NA)
  expect_equal(resolve_vector_ranks(vec), "A/B")
  expect_equal(resolve_vector_ranks(vec, strict = TRUE), "A/B/NA")
  expect_equal(
    resolve_vector_ranks(vec, method = "preference", preference_index = 1),
    "A"
  )
  expect_equal(
    resolve_vector_ranks(vec, method = "preference", preference_index = 2),
    "B"
  )

  expect_equal(
    resolve_vector_ranks(vec, method = "preference", preference_index = 3),
    "A/B"
  )
  expect_equal(
    resolve_vector_ranks(
      vec,
      method = "preference",
      preference_index = 3,
      second_method = c("rel_majority")
    ),
    "A/B"
  )

  expect_true(is.na(resolve_vector_ranks(
    vec,
    method = "preference",
    preference_index = 3,
    second_method = c("unanimity")
  )))
  expect_true(is.na(resolve_vector_ranks(
    vec,
    method = "preference",
    preference_index = 3,
    second_method = c("abs_majority")
  )))

  expect_error(resolve_vector_ranks(
    vec,
    method = "preference",
    preference_index = 4
  ))

  expect_true(is.na(resolve_vector_ranks(vec, method = "abs_majority")))

  expect_equal(resolve_vector_ranks(vec, method = "rel_majority"), "A/B")
  expect_true(is.na(resolve_vector_ranks(vec, method = "unanimity")))
  expect_true(is.na(resolve_vector_ranks(
    vec,
    method = "unanimity",
    strict = TRUE
  )))
})


test_that("resolve_vector_ranks works with a vector of one A, two B", {
  vec <- c("A", "B", "B")
  expect_equal(resolve_vector_ranks(vec), "A/B")
  expect_equal(
    resolve_vector_ranks(vec, method = "preference", preference_index = 1),
    "A"
  )
  expect_equal(resolve_vector_ranks(vec, method = "abs_majority"), "B")
  expect_equal(resolve_vector_ranks(vec, method = "rel_majority"), "B")

  expect_true(is.na(resolve_vector_ranks(
    vec,
    method = "abs_majority",
    nb_agree_threshold = 3
  )))
  expect_true(is.na(resolve_vector_ranks(
    vec,
    method = "rel_majority",
    nb_agree_threshold = 3
  )))
  expect_true(is.na(resolve_vector_ranks(vec, method = "unanimity")))
})

test_that("rename_samples works", {
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
  suppressWarnings(result <- is_mumu_installed())
  expect_type(result, "logical")
})

test_that("is_vsearch_installed works", {
  result <- is_vsearch_installed()
  expect_type(result, "logical")
})

test_that("is_cutadapt_installed works", {
  suppressWarnings(result <- is_cutadapt_installed())
  expect_type(result, "logical")
})

test_that("is_swarm_installed works", {
  result <- is_swarm_installed()
  expect_type(result, "logical")
})
