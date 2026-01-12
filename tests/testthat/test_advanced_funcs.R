data(data_fungi)

test_that("lulu works", {
  skip_on_cran()
  res1 <- lulu_pq(data_fungi_sp_known)
  expect_equal(length(res1), 4)
  expect_equal(ntaxa(res1$new_physeq), 549)
  expect_equal(nsamples(res1$new_physeq), 185)
  
  res2 <- lulu_pq(data_fungi_sp_known, verbose=TRUE, clean_pq=TRUE)
  expect_equal(length(res2), 4)
  expect_equal(ntaxa(res2$new_physeq), 549)
  expect_equal(nsamples(res2$new_physeq), 184)
})

test_that("glmutli_pq works", {
  if (requireNamespace("glmulti", quietly = TRUE)) {
    res_glmulti <- 
      glmutli_pq(data_fungi, "Hill_0 ~ Hill_1 + Abundance + Time + Height", level = 1)
    expect_equal(dim(res_glmulti), c(5, 6))
     res_glmulti_interaction <- 
       glmutli_pq(data_fungi, "Hill_0 ~ Abundance + Time + Height", level = 2)
    expect_equal(dim(res_glmulti_interaction), c(11, 6))
  }
})
