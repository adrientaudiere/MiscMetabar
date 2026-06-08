skip_on_cran()
data("data_fungi_mini", package = "MiscMetabar")

test_that("query_cov and target_cov args are exposed", {
  expect_true(all(
    c("query_cov", "target_cov") %in% names(formals(vsearch_clustering))
  ))
  expect_true(all(
    c("query_cov", "target_cov") %in% names(formals(postcluster_pq))
  ))
  expect_null(formals(vsearch_clustering)$query_cov)
  expect_null(formals(vsearch_clustering)$target_cov)
})

if (is_vsearch_installed()) {
  test_that("vsearch_clustering accepts query_cov and target_cov", {
    expect_s4_class(
      vsearch_clustering(
        data_fungi_mini,
        query_cov = 0.9,
        target_cov = 0.9
      ),
      "phyloseq"
    )
  })

  test_that("postcluster_pq passes query_cov and target_cov to vsearch", {
    expect_s4_class(
      postcluster_pq(
        data_fungi_mini,
        method = "vsearch",
        query_cov = 0.9,
        target_cov = 0.9
      ),
      "phyloseq"
    )
  })
}
