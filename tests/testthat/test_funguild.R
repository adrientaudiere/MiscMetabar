data(data_fungi)

withr::local_envvar(
  R_USER_CACHE_DIR = tempfile(),
  .local_envir = teardown_env()
)

test_that("funguild_assign works", {
  skip_on_cran()
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    # Get the database first
    db <- get_funguild_db()
    # Test with a simple taxonomy
    result <- funguild_assign(tax_table(data_fungi)[1:5, ], db)
    expect_type(result, "list")
  }
})

test_that("get_funguild_db works", {
  skip_on_cran()
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    db <- get_funguild_db()
    expect_s3_class(db, "data.frame")
  }
})

test_that("formattable_pq works", {
  skip_on_cran()
  if (requireNamespace("formattable", quietly = TRUE)) {
    result <- formattable_pq(data_fungi)
    expect_s3_class(result, "formattable")
  }
})
