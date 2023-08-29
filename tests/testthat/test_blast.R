# blast_to_phyloseq
# blast_pq
# filter_asv_blast
# add_blast_info
# blast_to_derep

data("data_fungi")
df_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
df_basidio <- subset_taxa_pq(df_basidio, colSums(df_basidio@otu_table) > 1000)

suppressWarnings(blast_error_or_not <- try(system("blastn 2>&1", intern = TRUE), silent = TRUE))

if (class(blast_error_or_not) == "try-error") {
  test_that("blast_pq send an error when Blast is not installed", {
    # expect_error()
  })
} else {
  test_that("blast_to_phyloseq works fine", {
    expect_s3_class(blast_on_df <- blast_to_phyloseq(df_basidio, "inst/extdata/1000_sp_UNITE_sh_general_release_dynamic.fasta"), "data.frame")
  })
  test_that("blast_pq works fine", {
    expect_s3_class(blast_df <- blast_pq(df_basidio, "inst/extdata/1000_sp_UNITE_sh_general_release_dynamic.fasta"), "data.frame")
  })

}
