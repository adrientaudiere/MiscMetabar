data(data_fungi)
data(enterotype)

test_that("filt_taxa_pq works", {
 
  result <- filt_taxa_pq(data_fungi, min_nb_seq = 20)
  expect_s4_class(result, "phyloseq")
  expect_equal(ntaxa(result), 1388)
 
  result <- filt_taxa_pq(data_fungi, min_occurence = 2)
    expect_s4_class(result, "phyloseq")
  expect_equal(ntaxa(result), 1214)

  result <- filt_taxa_pq(data_fungi,
  min_occurence = 2,
  min_nb_seq = 10, clean_pq = FALSE
)
  expect_s4_class(result, "phyloseq")
  expect_equal(ntaxa(result),1087)
  
result <- filt_taxa_pq(data_fungi,
  min_occurence = 4,
  min_nb_seq = 100,
  combination = "OR"
)
  expect_s4_class(result, "phyloseq")
  expect_equal(ntaxa(result), 1261)
})

test_that("filt_taxa_wo_NA works", {
  result <- filt_taxa_wo_NA(data_fungi, "Phylum")
  expect_s4_class(result, "phyloseq")
  expect_equal(ntaxa(result), 1327)
})

test_that("distri_1_taxa works", {
  result <- distri_1_taxa(data_fungi, "Height", "ASV2")
  expect_s3_class(result, "data.frame")
})
