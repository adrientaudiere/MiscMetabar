path_db <- system.file("extdata",
  "100_sp_UNITE_sh_general_release_dynamic.fasta",
  package = "MiscMetabar", mustWork = TRUE
)

suppressWarnings(blast_error_or_not <-
  try(system("blastn 2>&1", intern = TRUE), silent = TRUE))

if (inherits(blast_error_or_not, "try-error")) {
  message(
    "blast_to_phyloseq(), filter_asv_blast(), blast_to_derep(),
          add_blast_info, and blast_pq() can't be tested when
          vsearch is not installed"
  )
} else {
  test_that(
    "blast_to_phyloseq works fine",
    {
      expect_s3_class(
        blast_on_df <-
          blast_to_phyloseq(data_fungi_mini, path_db),
        "data.frame"
      )
      expect_s3_class(
        blast_to_phyloseq(data_fungi_mini, path_db, list_no_output_query = TRUE),
        "data.frame"
      )
      expect_s3_class(
        blast_to_phyloseq(
          data_fungi_mini,
          path_db,
          unique_per_seq = TRUE,
          score_filter = FALSE
        ),
        "data.frame"
      )
      expect_error(blast_to_phyloseq(data_fungi_mini, "inst/extdata/nil.fasta"))
    }
  )

  test_that(
    "blast_pq works fine",
    {
      expect_s3_class(blast_df <-
        blast_pq(data_fungi_mini, path_db), "data.frame")
      expect_equal(ncol(blast_df), 9)
      expect_true(nrow(blast_df) > 0)
      expect_s3_class(
        blast_df <- blast_pq(data_fungi_mini, path_db,
          unique_per_seq = TRUE
        ),
        "data.frame"
      )
      expect_s3_class(
        blast_df <- blast_pq(data_fungi_mini, path_db,
          score_filter = FALSE
        ),
        "data.frame"
      )
      expect_s3_class(
        blast_df <-
          blast_pq(
            data_fungi_mini,
            path_db,
            unique_per_seq = FALSE,
            score_filter = FALSE
          ),
        "data.frame"
      )
      expect_s3_class(
        blast_df <-
          blast_pq(data_fungi_mini, database = "inst/extdata/dbase"),
        "data.frame"
      )
      expect_error(
        blast_df <-
          blast_pq(data_fungi_mini, fasta_for_db = path_db, database = "inst/extdata/dbase")
      )
      expect_error(blast_df <- blast_pq(data_fungi_mini))
    }
  )

  test_that(
    "filter_asv_blast works fine",
    {
      expect_s4_class(
        df_blast <-
          filter_asv_blast(data_fungi_mini, path_db),
        "phyloseq"
      )
      expect_s4_class(
        df_blast <-
          filter_asv_blast(
            data_fungi_mini,
            path_db,
            id_filter = 50,
            e_value_filter = 10,
            bit_score_filter = 20,
            min_cover_filter = 20
          ),
        "phyloseq"
      )
      expect_s4_class(
        df_blast <-
          filter_asv_blast(
            data_fungi_mini,
            path_db,
            add_info_to_taxtable = FALSE,
            nproc = 2
          ),
        "phyloseq"
      )
      expect_error(filter_asv_blast(data_fungi_mini, "inst/extdata/nil.fasta"))
    }
  )

  test_that(
    "filter_asv_blast returns NULL with strict filters",
    {
      expect_null(
        expect_message(
          filter_asv_blast(
            data_fungi_mini,
            path_db,
            id_filter = 100,
            e_value_filter = 0,
            bit_score_filter = 10000,
            min_cover_filter = 100
          ),
          "No taxa passed the filter criteria"
        )
      )
    }
  )

  test_that(
    "add_blast_info works fine",
    {
      expect_s4_class(
        df_blast <-
          add_blast_info(data_fungi_mini, path_db),
        "phyloseq"
      )
      expect_error(add_blast_info(data_fungi_mini, "inst/extdata/nil.fasta"))
    }
  )

  derep_data <-
    derepFastq(unlist(list_fastq_files("inst/extdata/")))
  test_that(
    "blast_to_derep works fine",
    {
      expect_s3_class(
        derep_blast <-
          blast_to_derep(derep_data, path_db),
        "data.frame"
      )
      expect_s3_class(
        derep_blast <-
          blast_to_derep(derep_data, path_db, score_filter = TRUE),
        "data.frame"
      )
      expect_s3_class(
        derep_blast2 <-
          blast_to_derep(derep_data, "inst/extdata/ex_little.fasta"),
        "data.frame"
      )
      expect_s3_class(
        derep_blast2 <-
          blast_to_derep(
            derep_data,
            "inst/extdata/ex_little.fasta",
            unique_per_seq = TRUE,
            score_filter = FALSE
          ),
        "data.frame"
      )
      expect_s3_class(
        derep_blast2 <-
          blast_to_derep(
            derep_data,
            "inst/extdata/ex_little.fasta",
            list_no_output_query = TRUE
          ),
        "data.frame"
      )
      expect_error(blast_to_derep(derep_data, "inst/extdata/nil.fasta"))
    }
  )
}

file.remove(list.files("tests/testthat", pattern = "dbase"))
