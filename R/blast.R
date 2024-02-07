################################################################################
#' Blast some sequence against `refseq` slot of a \code{\link{phyloseq-class}}
#'   object.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object obtained
#'   using the `phyloseq` package.
#' @param seq2search (required) path to a fasta file defining the sequences
#'   you want to blast against the ASV sequences from the physeq object.
#' @param blastpath path to blast program
#' @param id_cut (default: 90) cut of in identity percent to keep result
#' @param bit_score_cut (default: 50) cut of in bit score to keep result
#'   The higher the bit-score, the better the sequence similarity.
#'   The bit-score is the requires size of a sequence database in which the current
#'   match could be found just by chance. The bit-score is a log2 scaled and
#'   normalized raw-score. Each increase by one doubles the required database size
#'   (2bit-score).
#' @param min_cover_cut (default: 50) cut of in query cover (%) to keep result
#' @param e_value_cut (default: 1e-30)  cut of in e-value (%) to keep result
#'   The BLAST E-value is the number of expected hits of similar quality (score)
#'   that could be found just by chance.
#' @param unique_per_seq (logical, default FALSE) if TRUE only return the better match
#'  (higher **bit score**) for each sequence
#' @param score_filter (logical, default TRUE) does results are filter by score? If
#'   FALSE, `id_cut`,`bit_score_cut`, `e_value_cut` and `min_cover_cut` are ignored
#' @param list_no_output_query (logical) does the result table include
#'   query sequences for which `blastn` does not find any correspondence?
#' @param args_makedb Additional parameters parse to makeblastdb command
#' @param args_blastn Additional parameters parse to blastn command
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for blast (args -num_threads
#'   for blastn command)
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - db.fasta (refseq transformed into a database)
#'   - dbase list of files (output of blastn)
#'   - blast_result.txt the summary result of blastn using
#'     `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`
#'
#' @seealso  [MiscMetabar::blast_pq()] to use `refseq` slot as query sequences
#'   against un custom database.
#'
#' @return  the blast table
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' blastpath <- "...YOUR_PATH_TO_BLAST..."
#' blast_to_phyloseq(data_fungi,
#'   seq2search = system.file("extdata", "ex.fasta",
#'     package = "MiscMetabar", mustWork = TRUE
#'   ),
#'   blastpath = blastpath
#' )
#' }
blast_to_phyloseq <- function(physeq,
                              seq2search,
                              blastpath = NULL,
                              id_cut = 90,
                              bit_score_cut = 50,
                              min_cover_cut = 50,
                              e_value_cut = 1e-30,
                              unique_per_seq = FALSE,
                              score_filter = TRUE,
                              list_no_output_query = FALSE,
                              args_makedb = NULL,
                              args_blastn = NULL,
                              nproc = 1,
                              keep_temporary_files = FALSE) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, paste0(tempdir(), "/", "db.fasta"))

  system(paste0(
    blastpath,
    "makeblastdb -dbtype nucl -in ",
    paste0(tempdir(), "/", "db.fasta"),
    " -out ",
    paste0(tempdir(), "/", "dbase"),
    " ",
    args_makedb
  ))

  system(
    paste0(
      blastpath,
      "blastn -query ",
      seq2search,
      " -db ",
      paste0(tempdir(), "/", "dbase"),
      " -out ",
      paste0(tempdir(), "/", "blast_result.txt"),
      " -outfmt \"6 qseqid qlen sseqid slen",
      " length pident evalue bitscore qcovs\"",
      " -num_threads ", nproc,
      " ",
      args_blastn
    )
  )
  if (file.info(paste0(tempdir(), "/", "blast_result.txt"))$size > 0) {
    blast_tab <- utils::read.table(
      paste0(tempdir(), "/", "blast_result.txt"), ,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    blast_tab_OK <- TRUE
  } else {
    blast_tab_OK <- FALSE
  }

  if (!keep_temporary_files) {
    unlink(paste0(tempdir(), "/", "blast_result.txt"))
    unlink(list.files(tempdir(), pattern = "dbase"))
    unlink(paste0(tempdir(), "/", "db.fasta"))
  } else {
    message(paste0("Temporary files are located at ", tempdir()))
  }

  if (!blast_tab_OK) {
    message("None query sequences matched your phyloseq references sequences.")
    return(NULL)
  }

  names(blast_tab) <- c(
    "Query name",
    "Query seq. length",
    "Taxa name",
    "Taxa seq. length",
    "Alignment length",
    "% id. match",
    "e-value",
    "bit score",
    "Query cover"
  )

  blast_tab <- blast_tab[order(blast_tab[, "bit score"], decreasing = TRUE), ]

  if (unique_per_seq) {
    blast_tab <- blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if (score_filter) {
    blast_tab <- blast_tab[blast_tab[, "bit score"] > bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
    blast_tab <- blast_tab[blast_tab[, "Query cover"] > min_cover_cut, ]
    blast_tab <- blast_tab[blast_tab[, "e-value"] < e_value_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  if (list_no_output_query) {
    fasta_file <- Biostrings::readDNAStringSet(seq2search)
    seq_name <- names(fasta_file)
    no_output_query <- seq_name[!seq_name %in% blast_tab[1, ]]
    if (length(no_output_query) > 0) {
      mat_no_output_query <- matrix(NA,
        ncol = ncol(blast_tab),
        nrow = length(no_output_query)
      )
      mat_no_output_query[, 1] <- no_output_query
      colnames(mat_no_output_query) <- colnames(blast_tab)
      blast_tab <- rbind(blast_tab, mat_no_output_query)
    }
  } else {
    blast_tab <- blast_tab
  }

  if (nrow(blast_tab) == 0) {
    stop("No blast query match the score filters")
  }
  return(blast_tab)
}
################################################################################


################################################################################
#' Blast all sequence of `refseq` slot of a \code{\link{phyloseq-class}}
#'   object against a custom database.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams blast_to_phyloseq
#' @param fasta_for_db path to a fasta file to make the blast database
#' @param database path to a blast database
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - db.fasta (refseq transformed into a database)
#'   - dbase list of files (output of blastn)
#'   - blast_result.txt the summary result of blastn using
#'     `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`
#'
#' @seealso  [MiscMetabar::blast_to_phyloseq()] to use `refseq`
#'   slot as a database
#' @return  a blast table
#' @export
#'
blast_pq <- function(physeq,
                     fasta_for_db = NULL,
                     database = NULL,
                     blastpath = NULL,
                     id_cut = 90,
                     bit_score_cut = 50,
                     min_cover_cut = 50,
                     e_value_cut = 1e-30,
                     unique_per_seq = FALSE,
                     score_filter = TRUE,
                     nproc = 1,
                     args_makedb = NULL,
                     args_blastn = NULL,
                     keep_temporary_files = FALSE) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(
    dna,
    paste0(tempdir(), "/", "physeq_refseq.fasta")
  )

  if (is.null(fasta_for_db) && is.null(database)) {
    stop("The function required a value for the parameters
         `fasta_for_db` or `database` to run.")
  } else if (!is.null(fasta_for_db) && !is.null(database)) {
    stop("You assign values for both `fasta_for_db` and
         `database` args. Please use only one.")
  } else if (!is.null(fasta_for_db) && is.null(database)) {
    message("Build the database from fasta_for_db")
    system(paste0(
      blastpath,
      "makeblastdb -dbtype nucl -in ",
      fasta_for_db,
      " -out ",
      paste0(tempdir(), "/", "dbase"),
      " ",
      args_makedb
    ))
    message("Blast refseq from physeq object against the database")
    system(
      paste0(
        blastpath,
        "blastn -query ",
        paste0(tempdir(), "/", "physeq_refseq.fasta"),
        " -db ",
        paste0(tempdir(), "/", "dbase"),
        " -out ",
        paste0(tempdir(), "/", "blast_result.txt"),
        " -num_threads ", nproc,
        " -outfmt \"6 qseqid qlen sseqid slen",
        " length pident evalue bitscore qcovs\"",
        " ",
        args_blastn
      )
    )
  } else if (is.null(fasta_for_db) && !is.null(database)) {
    message("Blast refseq from physeq object against the database")
    system(
      paste0(
        blastpath,
        "blastn -query ",
        paste0(tempdir(), "/", "physeq_refseq.fasta"),
        " -db ",
        database,
        " -out ",
        paste0(tempdir(), "/", "blast_result.txt"),
        " -outfmt \"6 qseqid qlen sseqid slen",
        " length pident evalue bitscore qcovs\"",
        " -num_threads ", nproc,
        " ",
        args_blastn
      )
    )
  }

  if (file.info(paste0(tempdir(), "/", "blast_result.txt"))$size > 0) {
    blast_tab <- utils::read.table(
      paste0(tempdir(), "/", "blast_result.txt"), ,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    blast_tab_OK <- TRUE
  } else {
    blast_tab_OK <- FALSE
  }

  if (!keep_temporary_files) {
    unlink(paste0(tempdir(), "/", "blast_result.txt"))
    unlink(list.files(tempdir(), pattern = "dbase"))
    unlink(paste0(tempdir(), "/", "db.fasta"))
  } else {
    message(paste0("Temporary files are located at ", tempdir()))
  }

  if (!blast_tab_OK) {
    message("None query sequences matched your phyloseq references sequences.")
    return(NULL)
  }

  names(blast_tab) <- c(
    "Query name",
    "Query seq. length",
    "Taxa name",
    "Taxa seq. length",
    "Alignment length",
    "% id. match",
    "e-value",
    "bit score",
    "Query cover"
  )

  blast_tab <- blast_tab[order(blast_tab[, "bit score"], decreasing = TRUE), ]

  if (unique_per_seq) {
    blast_tab <- blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if (score_filter) {
    blast_tab <- blast_tab[blast_tab[, "bit score"] > bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
    blast_tab <- blast_tab[blast_tab[, "Query cover"] > min_cover_cut, ]
    blast_tab <- blast_tab[blast_tab[, "e-value"] < e_value_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  if (nrow(blast_tab) == 0) {
    stop("No blast query match the score filters")
  }
  return(blast_tab)
}



################################################################################
#' Filter undesirable taxa using blast against a against a custom database.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams blast_to_phyloseq
#' @param fasta_for_db path to a fasta file to make the blast database
#' @param database path to a blast database
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples and empty ASV are discarded
#'   after filtering.
#' @param id_filter (default: 90) cut of in identity percent to keep result
#' @param bit_score_filter (default: 50) cut of in bit score to keep result
#'   The higher the bit-score, the better the sequence similarity.
#'   The bit-score is the requires size of a sequence database in which the current
#'   match could be found just by chance. The bit-score is a log2 scaled and
#'   normalized raw-score. Each increase by one doubles the required database size
#'   (2bit-score).
#' @param min_cover_filter (default: 50) cut of in query cover (%) to keep result
#' @param e_value_filter (default: 1e-30)  cut of in e-value (%) to keep result
#'   The BLAST E-value is the number of expected hits of similar quality (score)
#'   that could be found just by chance.
#' @param add_info_to_taxtable (logical, default TRUE) Does the blast information
#'   are added to the taxtable ?
#' @param ... Others options for the `blast_pq()` function. See `?blast_pq`.
#'   Note that params `unique_per_seq` must be lft to TRUE and `score_filter`
#'   must be left to FALSE.
#' @export
#' @return A new \code{\link{phyloseq-class}} object.


filter_asv_blast <- function(physeq,
                             fasta_for_db = NULL,
                             database = NULL,
                             clean_pq = TRUE,
                             add_info_to_taxtable = TRUE,
                             id_filter = 90,
                             bit_score_filter = 50,
                             min_cover_filter = 50,
                             e_value_filter = 1e-30,
                             ...) {
  blast_tab <- blast_pq(
    physeq = physeq,
    fasta_for_db = fasta_for_db,
    database = database,
    unique_per_seq = TRUE,
    score_filter = FALSE,
    ...
  )

  condition <- blast_tab[, "Query cover"] > min_cover_filter &
    blast_tab[, "bit score"] > bit_score_filter &
    blast_tab[, "% id. match"] > id_filter &
    blast_tab[, "e-value"] < e_value_filter

  names(condition) <- blast_tab[, "Query name"]

  new_physeq <- subset_taxa_pq(physeq, condition, clean_pq = FALSE)

  if (clean_pq) {
    new_physeq <- clean_pq(new_physeq)
  }

  if (add_info_to_taxtable) {
    info_to_taxtable <- blast_tab %>%
      group_by(`Query name`) %>%
      slice(which.min(`e-value`)) %>%
      ungroup()
    new_physeq@tax_table <- tax_table(as.matrix(cbind(
      new_physeq@tax_table,
      info_to_taxtable[
        match(
          taxa_names(new_physeq),
          info_to_taxtable[, "Query name"]$`Query name`
        ),
        c(
          "Query name", "Taxa name", "bit score",
          "% id. match", "Query cover", "e-value"
        )
      ]
    )))
  }

  return(new_physeq)
}


#' Blast some sequence against sequences from of a \code{\link{derep-class}}
#'   object.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams blast_to_phyloseq
#' @param derep The result of `dada2::derepFastq()`. A list of `derep-class`
#'   object.
#' @param seq2search (required) path to a fasta file defining the sequences
#'   you want to blast against the ASV sequences from the physeq object.
#' @param min_length_seq (default: 200) Removed sequences with less than
#'   `min_length_seq` from derep before blast. Set to 0 to discard filtering
#'    sequences by length.
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary
#'   files :
#'   - db.fasta (refseq transformed into a database)
#'   - dbase list of files (output of blastn)
#'   - blast_result.txt the summary result of blastn using
#'     `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`
#' @return A blast table
#'
#' @export
#' @seealso  [MiscMetabar::blast_pq()] to use `refseq` slot as query sequences
#'   against un custom database and [MiscMetabar::blast_to_phyloseq()]  to use
#'    `refseq` slot as a database
#' @author Adrien Taudière
blast_to_derep <- function(derep,
                           seq2search,
                           blastpath = NULL,
                           id_cut = 90,
                           bit_score_cut = 50,
                           min_cover_cut = 50,
                           e_value_cut = 1e-30,
                           unique_per_seq = FALSE,
                           score_filter = FALSE,
                           list_no_output_query = FALSE,
                           min_length_seq = 200,
                           args_makedb = NULL,
                           args_blastn = NULL,
                           nproc = 1,
                           keep_temporary_files = FALSE) {
  if (!inherits(derep[[1]], "derep")) {
    stop("derep must be an object of class derep-class")
  }

  derep_list <- lapply(derep, function(l) {
    res <- names(l$uniques)
    res <- res[tapply(res, seq_along(res), nchar) > min_length_seq]
    return(res)
  })

  derep_occurence <- lapply(derep, function(l) {
    res <- names(l$uniques)
    res_occur <- as.vector(l$uniques)
    res <- res_occur[tapply(res, seq_along(res), nchar) > min_length_seq]
    return(res)
  })

  dna <- unlist(Biostrings::DNAStringSetList(derep_list))
  names(dna) <- paste0(names(dna), "(", unlist(derep_occurence), "seqs)")
  Biostrings::writeXStringSet(dna, paste0(tempdir(), "/", "db.fasta"))

  system(paste0(
    blastpath,
    "makeblastdb -dbtype nucl -in ",
    paste0(tempdir(), "/", "db.fasta"),
    " -out ",
    paste0(tempdir(), "/", "dbase"),
    " ", args_makedb
  ))

  system(
    paste0(
      blastpath,
      "blastn -query ",
      seq2search,
      " -db ",
      paste0(tempdir(), "/", "dbase"),
      " -out ",
      paste0(tempdir(), "/", "blast_result.txt"),
      " -outfmt \"6 qseqid qlen sseqid slen",
      " length pident evalue bitscore qcovs\"",
      " -num_threads ", nproc,
      " ",
      args_blastn
    )
  )


  if (file.info(paste0(tempdir(), "/", "blast_result.txt"))$size > 0) {
    blast_tab <- utils::read.table(
      paste0(tempdir(), "/", "blast_result.txt"), ,
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    blast_tab_OK <- TRUE
  } else {
    blast_tab_OK <- FALSE
  }

  if (!keep_temporary_files) {
    unlink(paste0(tempdir(), "/", "blast_result.txt"))
    unlink(list.files(tempdir(), pattern = "dbase"))
    unlink(paste0(tempdir(), "/", "db.fasta"))
  } else {
    message(paste0("Temporary files are located at ", tempdir()))
  }

  if (!blast_tab_OK) {
    message("None query sequences matched your phyloseq references sequences.")
    return(NULL)
  }

  names(blast_tab) <- c(
    "Query name",
    "Query seq. length",
    "Sample name",
    "Taxa seq. length",
    "Alignment length",
    "% id. match",
    "e-value",
    "bit score",
    "Query cover"
  )

  blast_tab$occurence <- sub("seqs\\)", "",
    sub(".*\\(", "",
      blast_tab$`Sample name`,
      perl = TRUE
    ),
    perl = TRUE
  )

  blast_tab <- blast_tab[order(blast_tab[, "bit score"], decreasing = TRUE), ]

  if (unique_per_seq) {
    blast_tab <- blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if (score_filter) {
    blast_tab <- blast_tab[blast_tab[, "bit score"] > bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
    blast_tab <- blast_tab[blast_tab[, "Query cover"] > min_cover_cut, ]
    blast_tab <- blast_tab[blast_tab[, "e-value"] < e_value_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  if (list_no_output_query) {
    fasta_file <- Biostrings::readDNAStringSet(seq2search)
    seq_name <- names(fasta_file)
    no_output_query <- seq_name[!seq_name %in% blast_tab[1, ]]
    if (length(no_output_query) > 0) {
      mat_no_output_query <- matrix(NA,
        ncol = ncol(blast_tab),
        nrow = length(no_output_query)
      )
      mat_no_output_query[, 1] <- no_output_query
      colnames(mat_no_output_query) <- colnames(blast_tab)
      blast_tab <- rbind(blast_tab, mat_no_output_query)
    }
  } else {
    blast_tab <- blast_tab
  }

  if (nrow(blast_tab) == 0) {
    stop("No blast query match the score filters")
  }
  return(blast_tab)
}


################################################################################
#' Add information from [blast_pq()] to the `tax_table` slot of a *phyloseq* object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Basically a wrapper of [blast_pq()] with option `unique_per_seq = TRUE` and
#'  `score_filter = FALSE`.
#'
#' Add the information to the taxtable
#'
#' @inheritParams clean_pq
#' @param fasta_for_db path to a fasta file to make the blast database
#' @param silent (logical) If true, no message are printing.
#' @param ... Other arguments passed on to [blast_pq()] function.
#' @return a physeq object with more information in tax_table based on a
#'   blast on a given database
#'
#' @export
#'
#' @author Adrien Taudière

add_blast_info <- function(physeq, fasta_for_db, silent = FALSE, ...) {
  verify_pq(physeq)
  res_blast <- blast_pq(physeq,
    fasta_for_db = fasta_for_db,
    unique_per_seq = TRUE, score_filter = FALSE, ...
  )
  new_physeq <- physeq

  new_physeq@tax_table <- tax_table(cbind(
    new_physeq@tax_table,
    as.matrix(res_blast[match(
      taxa_names(new_physeq),
      res_blast$`Query name`
    ), ])
  ))

  verify_pq(new_physeq)
  if (!silent) {
    message(paste0(
      "Add ", ncol(new_physeq@tax_table) - ncol(physeq@tax_table),
      " columns to taxonomic table"
    ))
  }
  return(new_physeq)
}
################################################################################
