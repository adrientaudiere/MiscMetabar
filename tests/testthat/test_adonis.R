data(data_fungi)

test_that("adonis function works fine", {
    expect_s3_class(adonis_pq(data_fungi, "Tree_name"), "anova")
    expect_s3_class(adonis_pq(subset_samples(data_fungi, !is.na(data_fungi@sam_data$Time)), "Time*Tree_name"), "anova")
    expect_error(adonis_pq(data_fungi, "Time*Tree_name"))
})