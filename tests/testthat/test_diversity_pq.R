skip_on_cran()
skip_if_not_installed("ecodive")

pq_10 <- prune_samples(sample_names(data_fungi_mini)[1:10], data_fungi_mini)
pq_10 <- clean_pq(pq_10, silent = TRUE)

test_that("diversity_pq() computes beta-diversity metrics matching vegan", {
  d_bray <- diversity_pq(pq_10, metric = "bray")
  expect_s3_class(d_bray, "dist")
  expect_equal(attr(d_bray, "Size"), nsamples(pq_10))

  otu <- as(otu_table(taxa_as_columns(pq_10)), "matrix")
  mode(otu) <- "numeric"
  d_veg <- vegan::vegdist(otu, method = "bray")
  expect_equal(as.numeric(d_bray), as.numeric(d_veg), tolerance = 1e-12)

  d_ait <- diversity_pq(pq_10, metric = "aitchison", pseudocount = 1)
  d_veg_ait <- vegan::vegdist(otu, method = "aitchison", pseudocount = 1)
  expect_equal(as.numeric(d_ait), as.numeric(d_veg_ait), tolerance = 1e-10)
})

test_that("diversity_pq() supports robust_aitchison with and without OptSpace imputation", {
  d_rait_eco <- diversity_pq(pq_10, metric = "robust_aitchison", impute = FALSE)
  expect_s3_class(d_rait_eco, "dist")

  otu <- as(otu_table(taxa_as_columns(pq_10)), "matrix")
  mode(otu) <- "numeric"
  rclr_zero <- vegan::decostand(otu, method = "rclr", impute = FALSE)
  rclr_zero[is.na(rclr_zero)] <- 0
  expect_equal(
    as.numeric(d_rait_eco),
    as.numeric(stats::dist(rclr_zero)),
    tolerance = 1e-12
  )

  d_rait_imp <- diversity_pq(pq_10, metric = "robust.aitchison", impute = TRUE)
  expect_s3_class(d_rait_imp, "dist")
  expect_equal(
    as.numeric(d_rait_imp),
    as.numeric(vegan::vegdist(otu, method = "robust.aitchison")),
    tolerance = 1e-12
  )
})

test_that("diversity_pq() returns a named list for multiple beta metrics", {
  res_list <- diversity_pq(pq_10, metric = c("bray", "robust_aitchison"))
  expect_type(res_list, "list")
  expect_named(res_list, c("bray", "robust_aitchison"))
  expect_s3_class(res_list$bray, "dist")
  expect_s3_class(res_list$robust_aitchison, "dist")
})

test_that("diversity_pq() computes alpha-diversity tibble and augments phyloseq", {
  alpha_tbl <- diversity_pq(pq_10, metric = c("shannon", "simpson", "observed"))
  expect_s3_class(alpha_tbl, "tbl_df")
  expect_equal(nrow(alpha_tbl), nsamples(pq_10))
  expect_named(alpha_tbl, c("Sample", "shannon", "simpson", "observed"))

  pq_aug <- diversity_pq(
    pq_10,
    metric = c("shannon", "observed"),
    add_to_phyloseq = TRUE
  )
  expect_s4_class(pq_aug, "phyloseq")
  expect_true(all(c("shannon", "observed") %in% colnames(pq_aug@sam_data)))
})

test_that("diversity_pq() validates unknown or mixed metrics and missing tree", {
  expect_error(diversity_pq(pq_10, metric = "not_a_metric"), "Unknown")
  expect_error(
    diversity_pq(pq_10, metric = c("bray", "shannon")),
    "Cannot mix"
  )
  expect_error(
    diversity_pq(pq_10, metric = "bray", add_to_phyloseq = TRUE),
    "alpha-diversity"
  )
  expect_error(
    diversity_pq(pq_10, metric = "unweighted_unifrac"),
    "phy_tree"
  )
})
