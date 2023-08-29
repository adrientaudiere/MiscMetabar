# blast_to_phyloseq
# blast_pq
# filter_asv_blast
# add_blast_info
# blast_to_derep

data("data_fungi")
df_basidio <- subset_taxa(data_fungi, Phylum == "Basidiomycota")
df_basidio <- subset_taxa_pq(df_basidio, colSums(df_basidio@otu_table) > 1000)
path_db <- "inst/extdata/1000_sp_UNITE_sh_general_release_dynamic.fasta"

suppressWarnings(blast_error_or_not <- try(system("blastn 2>&1", intern = TRUE), silent = TRUE))

if (class(blast_error_or_not) == "try-error") {
  test_that("blast_pq send an error when Blast is not installed", {
    # expect_error()
  })
} else {

  test_that("blast_to_phyloseq works fine", {
    expect_s3_class(blast_on_df <- blast_to_phyloseq(df_basidio, path_db), "data.frame")
    expect_error(blast_to_phyloseq(df_basidio, "inst/extdata/nil.fasta"))
  })

  test_that("blast_pq works fine", {
    expect_s3_class(blast_df <- blast_pq(df_basidio, path_db), "data.frame")
    expect_s3_class(blast_df <- blast_pq(df_basidio, path_db, list_no_output_query = TRUE), "data.frame")
    expect_s3_class(blast_df <- blast_pq(df_basidio, path_db, unique_per_seq = TRUE), "data.frame")
    expect_s3_class(blast_df <- blast_pq(df_basidio, path_db, score_filter = FALSE), "data.frame")
    expect_error(blast_pq(df_basidio, "inst/extdata/nil.fasta"))
  })

  test_that("filter_asv_blast works fine", {
    expect_s4_class(df_blast <- filter_asv_blast(df_basidio, path_db), "phyloseq")
    expect_s4_class(df_blast <- filter_asv_blast(df_basidio, path_db, unique_per_seq = TRUE), "phyloseq")
    expect_s4_class(df_blast <- filter_asv_blast(df_basidio, path_db, unique_per_seq = TRUE), "phyloseq")

    expect_error(filter_asv_blast(df_basidio, "inst/extdata/nil.fasta"))
  })

  test_that("add_blast_info works fine", {
    expect_s4_class(df_blast <- add_blast_info(df_basidio, path_db), "phyloseq")
    expect_error(add_blast_info(df_basidio, "inst/extdata/nil.fasta"))
  })

  test_that("add_blast_info works fine", {
    expect_s4_class(derep_blast <- blast_to_derep(df_basidio, path_db), "phyloseq")
    expect_s4_class(derep_blast <- blast_to_derep(df_basidio, ex.fasta), "phyloseq")
    expect_error(blast_to_derep(df_basidio, "inst/extdata/nil.fasta"))
  })
}
