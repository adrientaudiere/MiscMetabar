################################################################################
#' Blast some sequence against `refseq` slot of a \code{\link[phyloseq]{phyloseq-class}}
#'   object.
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Use the blast software.
#'
#' @param physeq (required): a \code{\link[phyloseq]{phyloseq-class}} object obtained
#'   using the `phyloseq` package.
#' @param seq2search (required) path to a fasta file defining the sequences
#'   you want to blast against the taxa (ASV, OTU) sequences from the physeq object.
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

  system(
    paste0(
      blastpath,
      "makeblastdb -dbtype nucl -in ",
      paste0(tempdir(), "/", "db.fasta"),
      " -out ",
      paste0(tempdir(), "/", "dbase"),
      " ",
      args_makedb
    )
  )

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
      " -num_threads ",
      nproc,
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

  if (keep_temporary_files) {
    message(paste0("Temporary files are located at ", tempdir()))
  } else {
    unlink(paste0(tempdir(), "/", "blast_result.txt"))
    unlink(list.files(tempdir(), pattern = "dbase"))
    unlink(paste0(tempdir(), "/", "db.fasta"))
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
#' Blast all sequence of `refseq` slot of a \code{\link[phyloseq]{phyloseq-class}}
#'   object against a custom database.
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use the blast software.
#'
#' @inheritParams blast_to_phyloseq
#' @param fasta_for_db Either a DNAStringSet object or a path to a fasta
#'  file to make the blast database.
#' @param database path to a blast database
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - db.fasta (refseq transformed into a database)
#'   - dbase list of files (output of blastn)
#'   - blast_result.txt the summary result of blastn using
#'     `-outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore qcovs"`
#'   - temp.fasta if fasta_for_db was a DNAStringSet object.
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
  Biostrings::writeXStringSet(dna, paste0(tempdir(), "/", "physeq_refseq.fasta"))

  if (inherits(fasta_for_db, "DNAStringSet")) {
    Biostrings::writeXStringSet(fasta_for_db, paste0(tempdir(), "/", "temp.fasta"))
    fasta_for_db <- paste0(tempdir(), "/", "temp.fasta")
  }

  if (is.null(fasta_for_db) && is.null(database)) {
    stop(
      "The function required a value for the parameters
         `fasta_for_db` or `database` to run."
    )
  } else if (!is.null(fasta_for_db) && !is.null(database)) {
    stop("You assign values for both `fasta_for_db` and
         `database` args. Please use only one.")
  } else if (!is.null(fasta_for_db) && is.null(database)) {
    message("Build the database from fasta_for_db")
    system(
      paste0(
        blastpath,
        "makeblastdb -dbtype nucl -in ",
        fasta_for_db,
        " -out ",
        paste0(tempdir(), "/", "dbase"),
        " ",
        args_makedb
      )
    )
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
        " -num_threads ",
        nproc,
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
        " -num_threads ",
        nproc,
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

  if (keep_temporary_files) {
    message(paste0("Temporary files are located at ", tempdir()))
  } else {
    unlink(paste0(tempdir(), "/", "blast_result.txt"))
    unlink(list.files(tempdir(), pattern = "dbase"))
    unlink(paste0(tempdir(), "/", "db.fasta"))
    unlink(paste0(tempdir(), "/", "temp.fasta"))
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
#' Filter undesirable taxa using blast against a custom database.
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use the blast software.
#'
#' @inheritParams blast_to_phyloseq
#' @param fasta_for_db path to a fasta file to make the blast database
#' @param database path to a blast database
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples and empty taxa (ASV, OTU) are discarded
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
#'   Note that params `unique_per_seq` must be left to TRUE and `score_filter`
#'   must be left to FALSE.
#' @export
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object.


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


################################################################################
#' @rdname filter_asv_blast
#' @export
filter_taxa_blast <- filter_asv_blast
################################################################################


#' Blast some sequence against sequences from of a \code{\link[dada2]{derep-class}}
#'   object.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use the blast software.
#'
#' @inheritParams blast_to_phyloseq
#' @param derep The result of `dada2::derepFastq()`. A list of `derep-class`
#'   object.
#' @param seq2search (required) path to a fasta file defining the sequences
#'   you want to blast against the taxa (ASV, OTU) sequences from the physeq object.
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

  if (keep_temporary_files) {
    message(paste0("Temporary files are located at ", tempdir()))
  } else {
    unlink(paste0(tempdir(), "/", "blast_result.txt"))
    unlink(list.files(tempdir(), pattern = "dbase"))
    unlink(paste0(tempdir(), "/", "db.fasta"))
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Basically a wrapper of [blast_pq()] with option `unique_per_seq = TRUE` and
#'  `score_filter = FALSE`.
#'
#' Add the information to the taxtable
#'
#' @inheritParams clean_pq
#' @param fasta_for_db path to a fasta file to make the blast database
#' @param silent (logical) If true, no message are printing.
#' @param suffix (character) The suffix to name the new columns.
#'   Set the suffix to "" in order to remove any suffix.
#' @param ... Other arguments passed on to [blast_pq()] function.
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with more information in tax_table based on a
#'   blast on a given database
#'
#' @export
#'
#' @author Adrien Taudière

add_blast_info <- function(physeq, fasta_for_db, silent = FALSE, suffix = "blast_info", ...) {
  verify_pq(physeq)
  res_blast <- blast_pq(physeq,
    fasta_for_db = fasta_for_db,
    unique_per_seq = TRUE,
    score_filter = FALSE,
    ...
  )
  new_physeq <- physeq

  new_taxtab <- cbind(
    new_physeq@tax_table,
    as.matrix(res_blast[match(
      taxa_names(new_physeq),
      res_blast$`Query name`
    ), ])
  )

  # colnames(new_taxtab) <- c(colnames(new_physeq@tax_table), paste0(colnames(physeq), suffix))

  new_physeq@tax_table <- tax_table(new_taxtab)

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


################################################################################
#' Assign taxonomy using blastn algorithm and the blast software
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use the blast software.
#'
#'
#' @inheritParams clean_pq
#' @param ref_fasta Either a DNAStringSet object or a path to a fasta
#'  file to make the blast database. It must be in sintax format.
#'  See [assign_sintax()].
#' @param database path to a blast database. Only used if ref_fasta
#'   is not set.
#' @param blastpath path to blast program.
#' @param behavior Either "return_matrix" (default), or "add_to_phyloseq":
#'
#'  - "return_matrix" return a list of two matrix with taxonomic value in the
#'    first element of the list and bootstrap value in the second one.
#'
#'  - "add_to_phyloseq" return a phyloseq object with amended slot `@taxtable`.
#'    Only available if using physeq input and not seq2search input.
#'
#' @param method (One of "vote" or "top-hit"). If top-hit, only the
#'   better match is used to assign taxonomy. If vote, the algorithm
#'   takes all (or `nb_voting` if `nb_voting` is not null) select assignation
#'   and resolve the conflict using the function [resolve_vector_ranks()].
#' @param suffix (character) The suffix to name the new columns.
#'   If set to "" (the default), the taxo_rank algorithm is used
#'   without suffix.
#' @param min_id (default: 95) the identity percent to take into account
#'   a references taxa
#' @param min_bit_score (default: 50) the minimum bit score to take
#'   into account a references taxa
#' @param min_cover (default: 50) cut of in query cover (%) to keep result
#' @param min_e_value (default: 1e-30)  cut of in e-value (%) to keep result
#'   The BLAST E-value is the number of expected hits of similar quality (score)
#'   that could be found just by chance.
#' @param nb_voting (Int, default NULL). The number of taxa to keep before apply
#'   a vote to resolve conflict. If NULL all taxa passing the filters (min_id,
#'   min_bit_score, min_cover and min_e_value) are selected.
#' @param column_names A vector of names for taxonomic ranks. Must
#'   correspond to names in the ref_fasta files.
#' @param vote_algorithm the method to vote among "consensus", "rel_majority",
#'  "abs_majority" and "unanimity". See [resolve_vector_ranks()] for more details.
#' @param strict (Logical, default FALSE). See [resolve_vector_ranks()] for more details.
#' @param nb_agree_threshold  See [resolve_vector_ranks()] for more details.
#' @param preference_index  See [resolve_vector_ranks()] for more details.
#' @param collapse_string  See [resolve_vector_ranks()] for more details.
#' @param replace_collapsed_rank_by_NA (Logical, default TRUE) See [resolve_vector_ranks()] for more details.
#' @param simplify_taxo (logical default TRUE). Do we apply the
#'   function [simplify_taxo()] to the phyloseq object?
#' @param keep_blast_metrics (Logical, default FALSE). If TRUE, the blast metrics
#'   ("Query seq. length", "Taxa seq. length", "Alignment length",  "% id. match", "e-value",
#'   "bit score" and "Query cover") are stored in the tax_table.
#' @param ... Other arguments passed on to [blast_pq()]
#'
#' @returns
#' - If behavior == "return_matrix" :
#'    - If method = "top-hit" a matrix of taxonomic assignation
#'    - If method = "vote", a list of two matrix, the first is the
#'      raw taxonomic assignation (before vote). The second one is
#'      the taxonomic assignation in which conflicts are resolved
#'      using vote.
#'
#' - If behavior == "vote", return a new phyloseq object
#' @export
#' @author Adrien Taudière
#' @examples
#' \dontrun{
#' ref_fasta <- Biostrings::readDNAStringSet(system.file("extdata",
#'   "mini_UNITE_fungi.fasta.gz",
#'   package = "MiscMetabar", mustWork = TRUE
#' ))
#'
#' # assign_blastn(data_fungi_mini, ref_fasta = ref_fasta) # error because not
#' # enough sequences in db so none blast query passed the filters.
#' # So we used low score filter hereafter.
#'
#' mat <- assign_blastn(data_fungi_mini,
#'   ref_fasta = ref_fasta,
#'   method = "top-hit", min_id = 70, min_e_value = 1e-3, min_cover = 50,
#'   min_bit_score = 20
#' )
#' head(mat)
#'
#' assign_blastn(data_fungi_mini,
#'   ref_fasta = ref_fasta, method = "vote",
#'   vote_algorithm = "rel_majority", min_id = 90, min_cover = 50,
#'   behavior = "add_to_phyloseq"
#' )@tax_table
#'
#' assign_blastn(data_fungi_mini,
#'   ref_fasta = ref_fasta, method = "vote",
#'   vote_algorithm = "consensus", replace_collapsed_rank_by_NA = FALSE,
#'   min_id = 90, min_cover = 50, behavior = "add_to_phyloseq"
#' )@tax_table
#' }
assign_blastn <- function(physeq,
                          ref_fasta = NULL,
                          database = NULL,
                          blastpath = NULL,
                          behavior = "return_matrix",
                          method = c("vote", "top-hit"),
                          suffix = "_blastn",
                          min_id = 95,
                          min_bit_score = 50,
                          min_cover = 95,
                          min_e_value = 1e-30,
                          nb_voting = NULL,
                          column_names = c(
                            "Kingdom",
                            "Phylum",
                            "Class",
                            "Order",
                            "Family",
                            "Genus",
                            "Species"
                          ),
                          vote_algorithm = c(
                            "consensus",
                            "rel_majority",
                            "abs_majority",
                            "unanimity"
                          ),
                          strict = FALSE,
                          nb_agree_threshold = 1,
                          preference_index = NULL,
                          collapse_string = "/",
                          replace_collapsed_rank_by_NA = TRUE,
                          simplify_taxo = TRUE,
                          keep_blast_metrics = FALSE,
                          ...) {
  method <- match.arg(method)
  vote_algorithm <- match.arg(vote_algorithm)
  if (method == "vote") {
    blast_tab_raw <- blast_pq(
      physeq = physeq,
      fasta_for_db = ref_fasta,
      database = database,
      unique_per_seq = FALSE,
      id_cut = min_id,
      bit_score_cut = min_bit_score,
      min_cover_cut = min_cover,
      e_value_cut = min_e_value,
      ...
    )

    if (is.null(nb_voting)) {
      nb_voting <- max(table(blast_tab_raw$`Query name`))
    }

    blast_tab <- blast_tab_raw |>
      tidyr::separate(`Taxa name`,
        into = c(paste0("Taxa_name_db", suffix), "Classification"),
        sep = ";tax="
      ) |>
      tidyr::separate("Classification",
        into = paste0(column_names, suffix),
        sep = ","
      ) |>
      group_by(`Query name`) |>
      slice_head(n = nb_voting) |>
      summarise(
        across(
          c(
            paste0(column_names, suffix),
            paste0("Taxa_name_db", suffix)
          ),
          ~ resolve_vector_ranks(
            .x,
            method = vote_algorithm,
            strict = strict,
            nb_agree_threshold = nb_agree_threshold,
            collapse_string = collapse_string,
            replace_collapsed_rank_by_NA = replace_collapsed_rank_by_NA
          )
        ),
        across(
          c(
            "Query seq. length",
            "Taxa seq. length",
            "Alignment length",
            "% id. match",
            "e-value",
            "bit score",
            "Query cover"
          ),
          list(mean = mean, sd = sd),
          .names = paste0("{.fn}_", "{.col}", suffix)
        )
      ) |>
      rename("taxa_names" = `Query name`)

    if (!keep_blast_metrics) {
      blast_tab <- blast_tab |>
        select(c(
          taxa_names,
          paste0(column_names, suffix),
          paste0("Taxa_name_db", suffix)
        ))
    }
  } else if (method == "top-hit") {
    blast_tab_raw <- blast_pq(
      physeq = physeq,
      fasta_for_db = ref_fasta,
      database = database,
      unique_per_seq = TRUE,
      id_cut = min_id,
      bit_score_cut = min_bit_score,
      min_cover_cut = min_cover,
      e_value_cut = min_e_value,
      ...
    )

    blast_tab <- blast_tab_raw |>
      tidyr::separate(`Taxa name`,
        into = c(paste0("Taxa_name_db", suffix), "Classification"),
        sep = ";tax="
      ) |>
      tidyr::separate("Classification",
        into = paste0(column_names, suffix),
        sep = ","
      ) |>
      rename("taxa_names" = `Query name`)

    if (!keep_blast_metrics) {
      blast_tab <- blast_tab |>
        select(c(
          taxa_names,
          paste0(column_names, suffix),
          paste0("Taxa_name_db", suffix)
        ))
    }
  }
  if (behavior == "return_matrix") {
    if (method == "vote") {
      return(list(
        "raw_blast_table" = blast_tab_raw,
        "blast_table_per_query" = blast_tab
      ))
    } else if (method == "top-hit") {
      return("blast_table_per_query" = blast_tab)
    }
  } else if (behavior == "add_to_phyloseq") {
    tax_tab <- as.data.frame(as.matrix(physeq@tax_table))
    tax_tab$taxa_names <- taxa_names(physeq)

    new_physeq <- physeq

    new_tax_tab <- left_join(tax_tab, blast_tab, by = join_by(taxa_names)) |>
      dplyr::select(-taxa_names) |>
      as.matrix()

    new_physeq@tax_table <- tax_table(new_tax_tab)
    taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

    if (simplify_taxo) {
      new_physeq <- simplify_taxo(new_physeq)
    }

    verify_pq(new_physeq)

    return(new_physeq)
  } else {
    stop("Param behavior must take either 'return_matrix' or 'add_to_phyloseq' value")
  }
}
################################################################################

################################################################################
#' Resolve conflict in a vector of taxonomy values
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Internaly used in the function [assign_blastn()] with method="vote"
#' @param vec (required) A vector of (taxonomic) values
#' @param method One of "consensus", "rel_majority", "abs_majority",
#'   "preference" or "unanimity". See details.
#'
#' @param strict (logical, default FALSE). If TRUE, NA are considered as
#'   informative in resolving conflict (i.e. NA are taking into account in vote).
#'   See details for more informations.
#' @param second_method One of "consensus", "rel_majority", "abs_majority",
#'   or "unanimity". Only used if method = "preference". See details.
#' @param nb_agree_threshold (Int, default 1) The minimum number of times a
#'   value must arise to be selected using vote. If 2, we only kept
#'   taxonomic value present at least 2 times in the vector.
#' @param preference_index (Int. default NULL). Required if method="preference".
#'   Useless for other method. The preference index is the index of the value in
#'   vec for which we have a preference.
#' @param collapse_string (default '/'). The character to collapse taxonomic names
#'   when multiple assignment is done.
#' @param replace_collapsed_rank_by_NA (logical, default FALSE). If set to TRUE,
#'   all multiple assignments (all taxonomic rank including the 'collapse_string'
#'   parameter) are replaced by NA.
#'
#' @details
#'
#' - `unanimity`: Only keep taxonomic value when all methods are agree
#' 	 - By default, the value with NA are not taking into account (strict=FALSE)
#' 	 - If `strict` , one NA in the row is sufficient to return a NA
#'
#' - `consensus`: Keep all taxonomic values separated by a '/' (separation can be modify using param `collapse_string`)
#' 	 - If `strict` is TRUE, NA are also added to taxonomic vector such as 'Tiger/Cat/NA' instead of 'Tiger/Cat'
#'
#' - `abs_majority`: Keep the most found taxonomic value if it represent at least half of all taxonomic values.
#' 	 - If `strict` is TRUE, NA values are also count to determine the majority. For example, a vector of taxonomic rank c("A", "A", "A", "B", NA, NA) will give a value of 'A' if `strict` is FALSE (default) but a value of NA if `strict` is TRUE.
#'
#' - `rel_majority`: Keep the most found taxonomic value. In case of equality, apply a consensus strategy (collapse values separated by a '/') across the most found taxonomic values.
#' 	 - If `strict` is TRUE, NA are considered as a rank and can win the relative majority vote.  If `strict` is FALSE (default), NA are removed before ranking the taxonomic values.
#' 	 - `nb_agree_threshold`: Only keep return value when at least x methods agreed with x is set by parameter `nb_agree_threshold`. By default, (`nb_agree_threshold` = 1): a majority of one is enough.
#'
#' - `preference`: Keep the value from a preferred column.
#' 	 - when the value is NA in the preferred column, apply a second strategy (by default `consensus`) to the other column (parameter `second_method`). Note that the parameters `strict` and `nb_agree_threshold` are used for the second_method consensus.
#' @returns a vector of length 1 (one character value)
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#'
#' resolve_vector_ranks(c("A"))
#' resolve_vector_ranks(c("A"),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A"), method = "abs_majority")
#' resolve_vector_ranks(c("A"), method = "rel_majority")
#' resolve_vector_ranks(c("A"),
#'   method = "rel_majority",
#'   nb_agree_threshold = 2
#' )
#' resolve_vector_ranks(c("A"), method = "unanimity")
#'
#' resolve_vector_ranks(c("A", "A", "A"))
#' resolve_vector_ranks(c("A", "A", "A"),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "A", "A"), method = "abs_majority")
#' resolve_vector_ranks(c("A", "A", "A"), method = "rel_majority")
#' resolve_vector_ranks(c("A", "A", "A"), method = "unanimity")
#'
#' resolve_vector_ranks(c(NA, NA, NA))
#' resolve_vector_ranks(c(NA, NA, NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c(NA, NA, NA), method = "abs_majority")
#' resolve_vector_ranks(c(NA, NA, NA), method = "rel_majority")
#' resolve_vector_ranks(c(NA, NA, NA), method = "unanimity")
#'
#' resolve_vector_ranks(c("A", "A", NA))
#' resolve_vector_ranks(c("A", "A", NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "A", NA), method = "rel_majority")
#' resolve_vector_ranks(c("A", "A", NA), method = "abs_majority")
#' resolve_vector_ranks(c("A", "A", NA, NA),
#'   method = "abs_majority",
#'   strict = FALSE
#' )
#' resolve_vector_ranks(c("A", "A", NA, NA),
#'   method = "abs_majority",
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "A", NA), method = "unanimity")
#' resolve_vector_ranks(c("A", "A", NA),
#'   method = "unanimity",
#'   strict = TRUE
#' )
#'
#' resolve_vector_ranks(c("A", "B", NA))
#' resolve_vector_ranks(c("A", "B", NA), strict = TRUE)
#' resolve_vector_ranks(c("A", "B", NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "B", NA), method = "abs_majority")
#' resolve_vector_ranks(c("A", "B", NA), method = "rel_majority")
#' resolve_vector_ranks(c("A", "B", NA),
#'   method = "rel_majority",
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "B", NA),
#'   method = "rel_majority",
#'   nb_agree_threshold = 2
#' )
#' resolve_vector_ranks(c("A", "B", NA), method = "unanimity")
#'
#' resolve_vector_ranks(c("A", NA, NA))
#' resolve_vector_ranks(c("A", NA, NA), method = "rel_majority")
#' resolve_vector_ranks(c("A", NA, NA), method = "unanimity")
#' resolve_vector_ranks(c("A", NA, NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", NA, NA),
#'   method = "preference",
#'   preference_index = 2
#' )
#' resolve_vector_ranks(c("A", NA, "B"),
#'   method = "preference",
#'   preference_index = 2
#' )
#' resolve_vector_ranks(c("A", NA, "B"),
#'   method = "preference",
#'   preference_index = 2, second_method = "abs_majority"
#' )
#'
#' resolve_vector_ranks(c("A", "B", "B"))
#' resolve_vector_ranks(c("A", "B", "B"),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "B", "B"), method = "abs_majority")
#' resolve_vector_ranks(c("A", "B", "B"), method = "rel_majority")
#' resolve_vector_ranks(c("A", "B", "B"), method = "unanimity")
#'
#'
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA))
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   method = "abs_majority"
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   method = "abs_majority",
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   method = "preference", preference_index = 6, second_method = "abs_majority"
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA, NA),
#'   method = "preference", preference_index = 6, second_method = "abs_majority"
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA, NA),
#'   method = "preference", preference_index = 6, second_method = "abs_majority",
#'   strict = TRUE
#' )
resolve_vector_ranks <- function(vec,
                                 method = c(
                                   "consensus",
                                   "rel_majority",
                                   "abs_majority",
                                   "preference",
                                   "unanimity"
                                 ),
                                 strict = FALSE,
                                 second_method = c(
                                   "consensus",
                                   "rel_majority",
                                   "abs_majority",
                                   "unanimity"
                                 ),
                                 nb_agree_threshold = 1,
                                 preference_index = NULL,
                                 collapse_string = "/",
                                 replace_collapsed_rank_by_NA = FALSE) {
  method <- match.arg(method)
  second_method <- match.arg(second_method)

  if (sum(is.na(vec)) == length(vec)) {
    res <- NA
    return(res)
  }

  if (method == "consensus") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    res <- paste0(unique(vec), collapse = collapse_string)
  } else if (method == "preference") {
    if (is.null(preference_index)) {
      stop("You must specify a preference_index if method = 'preference'")
    }
    res <- vec[preference_index]
    if (is.na(res)) {
      res <- resolve_vector_ranks(
        vec[-preference_index],
        method = second_method,
        nb_agree_threshold = nb_agree_threshold,
        strict = strict,
        collapse_string = collapse_string,
        replace_collapsed_rank_by_NA = replace_collapsed_rank_by_NA
      )
    }
  } else if (method == "abs_majority") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    nval <- length(vec)
    higher_val <- sort(table(vec), decreasing = T)[1]
    if (higher_val / nval > 0.5) {
      res <- names(higher_val)
    } else {
      res <- NA
    }
  } else if (method == "rel_majority") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    nval <- sum(table(vec, useNA = "ifany") == max(table(vec, useNA = "ifany")))
    if (sum(sort(table(vec, useNA = "ifany"), decreasing = T)[1:nval]) >= nb_agree_threshold) {
      res <- paste0(names(sort(table(vec, useNA = "ifany"), decreasing = T)[1:nval]), collapse = collapse_string)
    } else {
      res <- NA
    }
  } else if (method == "unanimity") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    if (length(unique(vec)) == 1) {
      res <- unique(vec)
    } else {
      res <- NA
    }
  }
  if (replace_collapsed_rank_by_NA &&
    sum(grepl(collapse_string, res)) > 0) {
    res <- NA
  }

  return(res)
}
################################################################################
