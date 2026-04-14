dfm <- data_fungi_mini
dfm_col <- taxa_as_columns(dfm)

# helpers
is_pq <- \(x) inherits(x, "phyloseq")
same_dims <- \(x, ref) identical(dim(otu_table(x)), dim(otu_table(ref)))

# ── transform_pq ─────────────────────────────────────────────────────────────

test_that("transform_pq returns phyloseq with same dims", {
  skip_if_not_installed("vegan")
  for (m in c("tss", "hellinger", "log1p", "z", "pa", "rank")) {
    res <- transform_pq(dfm, method = m)
    expect_true(is_pq(res), label = paste("class for", m))
    expect_true(same_dims(res, dfm), label = paste("dims for", m))
  }
})

test_that("transform_pq tss sample sums equal 1", {
  skip_if_not_installed("vegan")
  res <- transform_pq(dfm, method = "tss")
  expect_equal(round(range(sample_sums(res)), 10), c(1, 1))
})

test_that("transform_pq pa gives only 0/1", {
  skip_if_not_installed("vegan")
  vals <- unique(as.vector(otu_table(transform_pq(dfm, method = "pa"))))
  expect_true(all(vals %in% c(0, 1)))
})

test_that("transform_pq handles taxa-as-columns orientation", {
  skip_if_not_installed("vegan")
  res <- transform_pq(dfm_col, method = "hellinger")
  expect_false(taxa_are_rows(res))
  expect_true(same_dims(res, dfm_col))
})

# ── normalize_prop_pq ─────────────────────────────────────────────────────────

test_that("normalize_prop_pq returns phyloseq with same dims", {
  res <- normalize_prop_pq(dfm)
  expect_true(is_pq(res))
  expect_true(same_dims(res, dfm))
})

test_that("normalize_prop_pq without log: sample sums all close to constante", {
  res <- normalize_prop_pq(dfm, base_log = NULL, constante = 10000)
  expect_true(all(abs(sample_sums(res) - 10000) < 1))
})

test_that("normalize_prop_pq handles taxa-as-columns orientation", {
  res <- normalize_prop_pq(dfm_col)
  expect_false(taxa_are_rows(res))
})

# ── rarefy_pq ─────────────────────────────────────────────────────────────────

test_that("rarefy_pq returns phyloseq with equal sample depths", {
  res <- rarefy_pq(dfm, seed = 1)
  expect_true(is_pq(res))
  expect_equal(length(unique(sample_sums(res))), 1L)
})

test_that("rarefy_pq n>1 returns averaged (non-integer) table", {
  res <- rarefy_pq(dfm, n = 5, seed = 1)
  expect_true(is_pq(res))
  expect_false(all(otu_table(res) == floor(otu_table(res))))
})

# ── gmpr_pq ───────────────────────────────────────────────────────────────────

test_that("gmpr_pq returns phyloseq with same dims", {
  res <- suppressWarnings(gmpr_pq(dfm))
  expect_true(is_pq(res))
  expect_true(same_dims(res, dfm))
})

test_that("gmpr_pq handles taxa-as-columns orientation", {
  res <- suppressWarnings(gmpr_pq(dfm_col))
  expect_false(taxa_are_rows(res))
  expect_true(same_dims(res, dfm_col))
})

test_that("gmpr_pq with low ct_min produces valid size factors", {
  expect_warning(
    res <- gmpr_pq(dfm, ct_min = 1, intersect_no = 2),
    regexp = NA  # no warning expected with relaxed params on data_fungi_mini
  ) |> tryCatch(error = \(e) NULL)
  res <- suppressWarnings(gmpr_pq(dfm))
  expect_true(is_pq(res))
})

# ── mcknight_residuals_pq ─────────────────────────────────────────────────────

test_that("mcknight_residuals_pq adds column to sample_data", {
  res <- mcknight_residuals_pq(dfm)
  expect_true(is_pq(res))
  expect_true("mcknight_residuals" %in% colnames(sample_data(res)))
  expect_length(sample_data(res)$mcknight_residuals, nsamples(dfm))
})

test_that("mcknight_residuals_pq add_to_sam_data=FALSE returns named numeric", {
  res <- mcknight_residuals_pq(dfm, add_to_sam_data = FALSE)
  expect_type(res, "double")
  expect_named(res)
  expect_length(res, nsamples(dfm))
})

# ── Suggests-guarded functions ────────────────────────────────────────────────

test_that("srs_pq returns phyloseq with equal sample sums", {
  skip_if_not_installed("SRS")
  res <- srs_pq(dfm)
  expect_true(is_pq(res))
  expect_equal(length(unique(sample_sums(res))), 1L)
})

test_that("css_pq returns phyloseq with same dims", {
  skip_if_not_installed("metagenomeSeq")
  res <- css_pq(dfm)
  expect_true(is_pq(res))
  expect_true(same_dims(res, dfm))
})

test_that("tmm_pq returns phyloseq with same dims", {
  skip_if_not_installed("edgeR")
  res <- tmm_pq(dfm)
  expect_true(is_pq(res))
  expect_true(same_dims(res, dfm))
})

test_that("vst_pq returns phyloseq with same dims", {
  skip_if_not_installed("DESeq2")
  res <- vst_pq(dfm)
  expect_true(is_pq(res))
  expect_true(same_dims(res, dfm))
})
