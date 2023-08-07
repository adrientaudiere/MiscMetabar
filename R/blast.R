################################################################################
#' Blast some sequence against `refseq` slot of a \code{\link{phyloseq-class}}
#'   object.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @inheritParams clean_pq
#' @param seq2search (required) path to a fasta file defining the sequences
#'   you want to blast against the ASV sequences from the physeq object.
#' @param blastpath path to blast program
#' @param id_cut (default: 90) cut of in identity percent to keep result
#' @param bit_score_cut (default: 50) cut of in bit score to keep result
#' @param min_cover_cut (default: 50) cut of in query cover (%) to keep result
#' @param unique_per_seq (logical) if TRUE only return the first match for
#'   each sequence in seq2search
#' @param score_filter (logical) does results are filter by score? If
#'   FALSE, `id_cut`,`bit_score_cut` and `min_cover_cut` are ignored
#' @param list_no_output_query (logical) does the result table include
#'   query sequences for which `blastn` does not find any correspondence?
#' @param args_makedb Additional parameters parse to makeblastdb command
#' @param args_blastn Additional parameters parse to blastn command

#' @seealso  [MiscMetabar::blast_pq()] to use `refseq` slot as query sequences
#'   against un custom database.
#'
#' @return  the blast table
#' @export
#'
#' @examples
#' \dontrun{
#' data(data_fungi)
#' blastpath <- "...YOUR_PATH_TO_BLAST..."
#' blast_to_phyloseq(data_fungi,
#'   seq2search = "inst/extdata/ex.fasta",
#'   blastpath = blastpath
#' )
#' }
blast_to_phyloseq <- function(physeq,
                              seq2search,
                              blastpath = NULL,
                              id_cut = 90,
                              bit_score_cut = 50,
                              unique_per_seq = FALSE,
                              score_filter = TRUE,
                              list_no_output_query = FALSE,
                              args_makedb = NULL,
                              args_blastn = NULL) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, "db.fasta")

  system(paste(blastpath,
    "makeblastdb -dbtype nucl -in db.fasta -out dbase",
    " ", args_makedb,
    sep = ""
  ))

  system(
    paste(
      blastpath,
      "blastn -query ",
      seq2search,
      " -db dbase",
      " -out blast_result.txt",
      " -outfmt \"6 qseqid qlen sseqid slen",
      " length pident evalue bitscore qcovs\"",
      sep = ""
    )
  )

  if (file.info("blast_result.txt")$size > 0) {
    blast_tab <- utils::read.table(
      "blast_result.txt",
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    file.remove("blast_result.txt")
    file.remove(list.files(pattern = "dbase"))
    file.remove("db.fasta")
  } else {
    file.remove("blast_result.txt")
    file.remove(list.files(pattern = "dbase"))
    file.remove("db.fasta")
    stop("None query sequences matched your phyloseq references sequences.")
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

  blast_tab <- blast_tab[order(blast_tab[, "bit score"], decreasing = FALSE), ]

  if (unique_per_seq) {
    blast_tab <- blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if (score_filter) {
    blast_tab <- blast_tab[blast_tab[, "bit score"] > bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
    blast_tab <- blast_tab[blast_tab[, "Query cover"] > min_cover_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  if (list_no_output_query) {
    fastaFile <- Biostrings::readDNAStringSet(seq2search)
    seq_name <- names(fastaFile)
    no_output_query <- seq_name[!seq_name %in% blast_tab[1, ]]
    if (length(no_output_query) > 0) {
      mat_no_output_query <- matrix(NA,
        ncol = ncol(blast_tab),
        ncol = length(no_output_query)
      )
      mat_no_output_query[, 1] <- no_output_query
      colnames(mat_no_output_query) <- colnames(blast_tab)
      blast_tab <- rbind(blast_tab, mat_no_output_query)
    }
  } else {
    blast_tab <- blast_tab
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
#' @inheritParams clean_pq
#' @param fasta_for_db path to a fasta file to make the blast database
#' @param database path to a blast database
#' @param id_cut (default: 90) cut of in identity percent to keep result
#' @param bit_score_cut (default: 50) cut of in bit score to keep result
#' @param min_cover_cut (default: 50) cut of in query cover (%) to keep result
#' @param unique_per_seq (logical) if TRUE only return the first match for
#'   each sequence in seq2search
#' @param score_filter (logical) does results are filter by score? If
#'   FALSE, `id_cut`,`bit_score_cut` and `min_cover_cut` are ignored
#'  @param args_makedb Additional parameters parse to makeblastdb command
#'  @param args_blastn Additional parameters parse to blastn command
#'
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for blast (args -num_threads
#'   for blastn command)
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
                     unique_per_seq = FALSE,
                     score_filter = TRUE,
                     nproc = 1,
                     args_makedb = NULL,
                     args_blastn = NULL) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, "physeq_refseq.fasta")

  if (is.null(fasta_for_db) && is.null(database)) {
    stop("The function required a value for the parameters
         `fasta_for_db` or `database` to run.")
  } else if (!is.null(fasta_for_db) && !is.null(database)) {
    stop("You assign value for both `fasta_for_db` and
         `database` args. Please use only one.")
  } else if (!is.null(fasta_for_db) && is.null(database)) {
    message("Build the database from fasta_for_db")
    system(paste(blastpath,
      "makeblastdb -dbtype nucl -in ", fasta_for_db, " -out dbase",
      " ", args_makedb,
      sep = ""
    ))
    message("Blast refseq from physeq object against the database")
    system(
      paste(
        blastpath,
        "blastn -query ",
        "physeq_refseq.fasta",
        " -db dbase",
        " -out blast_result.txt",
        " -num_threads ", nproc,
        " -outfmt \"6 qseqid qlen sseqid slen",
        " length pident evalue bitscore qcovs\"",
        sep = ""
      )
    )
    if (file.info("blast_result.txt")$size > 0) {
      blast_tab <- utils::read.table(
        "blast_result.txt",
        sep = "\t",
        header = FALSE,
        stringsAsFactors = FALSE,
        comment.char = ""
      )
      file.remove("blast_result.txt")
      file.remove(list.files(pattern = "dbase"))
    } else {
      file.remove("blast_result.txt")
      file.remove(list.files(pattern = "dbase"))
      stop("None query sequences matched your phyloseq references sequences.")
    }
  } else if (is.null(fasta_for_db) && !is.null(database)) {
    message("Blast refseq from physeq object against the database")
    system(
      paste(
        blastpath,
        "blastn -query ",
        "physeq_refseq.fasta",
        " -db ", database,
        " -out blast_result.txt",
        " -num_threads ", nproc,
        " -outfmt \"6 qseqid qlen sseqid slen",
        " length pident evalue bitscore qcovs\"",
        sep = ""
      )
    )
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

  blast_tab <- blast_tab[order(blast_tab[, "% id. match"],
    decreasing = FALSE
  ), ]

  if (unique_per_seq) {
    blast_tab <- blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if (score_filter) {
    blast_tab <- blast_tab[blast_tab[, "bit score"] > bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
    blast_tab <- blast_tab[blast_tab[, "Query cover"] > min_cover_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  return(blast_tab)
}



################################################################################
#' Filter indesirable taxa using blast against a against a custom database.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams clean_pq
#' @param fasta_for_db path to a fasta file to make the blast database
#' @param database path to a blast database
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples and empty ASV are discarded
#'   after filtering.
#' @param blastpath path to blast program
#' @param id_cut (default: 80) cut of in identity percent to keep ASV
#' @param bit_score_cut (default: 150) cut of in bit score to keep result
#' @param min_cover_cut (default: 50) cut of in query cover (%) to keep result
#' @param add_info_to_taxtable: add some info from blast query to taxtable of the
#'   new physeq object. Only the information ("Query name", "Taxa name", "bit score",
#'   "% id. match", "Query cover", "e-value") for higher e-value hit.
#'   for each ASV is add to taxtable. Note that query name may be different from
#'   final taxa names as some function proposed to change the ASV names.
#' @param ... Others options for the `balst_pq()` function. See `?blast_pq`.
#'   Note that params `unique_per_seq` and `score_filter` must be set to TRUE.
#' @export
#' @return A new \code{\link{phyloseq-class}} object.


filter_asv_blast <- function(physeq,
                             fasta_for_db = NULL,
                             database = NULL,
                             clean_pq = TRUE,
                             blastpath = NULL,
                             id_cut = 80,
                             bit_score_cut = 150,
                             min_cover_cut = 50,
                             add_info_to_taxtable = TRUE,
                             ...) {
  blast_tab <- blast_pq(
    physeq = physeq,
    fasta_for_db = fasta_for_db,
    database = database,
    blastpath = blastpath,
    id_cut = id_cut,
    bit_score_cut = bit_score_cut,
    min_cover_cut = min_cover_cut,
    unique_per_seq = TRUE,
    score_filter = TRUE,
    ...
  )

  condition <- blast_tab[, "Query cover"] > min_cover_cut & blast_tab[, "bit score"] > bit_score_cut & blast_tab[, "% id. match"] > id_cut
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
        c("Query name", "Taxa name", "bit score", "% id. match", "Query cover", "e-value")
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
#' @param derep The result of `dada2::derepFastq()`. A list of `derep-class` object.
#' @param seq2search (required) path to a fasta file defining the sequences
#'   you want to blast against the ASV sequences from the physeq object.
#' @param blastpath path to blast program
#' @param id_cut (default: 90) cut of in identity percent to keep result
#' @param bit_score_cut (default: 50) cut of in bit score to keep result
#' @param min_cover_cut (default: 50) cut of in query cover (%) to keep result
#' @param unique_per_seq (logical) if TRUE only return the first match for
#'   each sequence in seq2search
#' @param score_filter (logical; default FALSE) does results are filter by score? If
#'   FALSE, `id_cut`,`bit_score_cut` and `min_cover_cut` are ignored.
#' @param list_no_output_query (logical) does the result table include
#'   query sequences for which `blastn` does not find any correspondence?
#' @param min_length_seq (default: 200) Removed sequences with less than
#'   `min_length_seq` from derep before blast. Set to 0 to discard filtering
#'    sequences by length.
#' @param args_makedb Additional parameters parse to makeblastdb command
#' @param args_blastn Additional parameters parse to blastn command
#'
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
                           unique_per_seq = FALSE,
                           score_filter = FALSE,
                           list_no_output_query = FALSE,
                           min_length_seq = 200,
                           args_makedb = NULL,
                           args_blastn = NULL) {
  if (!inherits(derep[[1]], "derep")) {
    stop("derep must be an object of class derep-class")
  }

  derep_list <- lapply(derep, function(l) {
    res <- names(l$uniques)
    res <- res[tapply(res, 1:length(res), nchar) > min_length_seq]
    return(res)
  })

  derep_occurence <- lapply(derep, function(l) {
    res <- names(l$uniques)
    res_occur <- as.vector(l$uniques)
    res <- res_occur[tapply(res, 1:length(res), nchar) > min_length_seq]
    return(res)
  })

  dna <- unlist(Biostrings::DNAStringSetList(derep_list))
  names(dna) <- paste0(names(dna), "(", unlist(derep_occurence), "seqs)")
  Biostrings::writeXStringSet(dna, "db.fasta")

  system(paste(
    blastpath,
    "makeblastdb -dbtype nucl -in db.fasta -out dbase",
    " ", args_makedb,
    sep = ""
  ))

  system(
    paste(
      blastpath,
      "blastn -query ",
      seq2search,
      " -db dbase",
      " -out blast_result.txt",
      " -outfmt \"6 qseqid qlen sseqid slen",
      " length pident evalue bitscore qcovs\"",
      sep = ""
    )
  )

  if (file.info("blast_result.txt")$size > 0) {
    blast_tab <- utils::read.table(
      "blast_result.txt",
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
    file.remove("blast_result.txt")
    file.remove(list.files(pattern = "dbase"))
    file.remove("db.fasta")
  } else {
    file.remove("blast_result.txt")
    file.remove(list.files(pattern = "dbase"))
    file.remove("db.fasta")
    stop("None query sequences matched your phyloseq references sequences.")
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

  blast_tab$occurence <- sub("seqs\\)", "", sub(".*\\(", "", blast_tab$`Sample name`, perl = T), perl = T)

  blast_tab <-
    blast_tab[order(blast_tab[, "e-value"], decreasing = FALSE), ]

  if (unique_per_seq) {
    blast_tab <- blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if (score_filter) {
    blast_tab <- blast_tab[blast_tab[, "bit score"] > bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
    blast_tab <-
      blast_tab[blast_tab[, "Query cover"] > min_cover_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  if (list_no_output_query) {
    fastaFile <- Biostrings::readDNAStringSet(seq2search)
    seq_name <- names(fastaFile)
    no_output_query <- seq_name[!seq_name %in% blast_tab[1, ]]
    if (length(no_output_query) > 0) {
      mat_no_output_query <- matrix(NA,
        ncol = ncol(blast_tab),
        ncol = length(no_output_query)
      )
      mat_no_output_query[, 1] <- no_output_query
      colnames(mat_no_output_query) <- colnames(blast_tab)
      blast_tab <- rbind(blast_tab, mat_no_output_query)
    }
  } else {
    blast_tab <- blast_tab
  }

  return(blast_tab)
}
