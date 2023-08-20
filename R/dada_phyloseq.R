if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

################################################################################
#' Add dna in `refseq` slot of a `physeq` object using taxa names and renames taxa
#'   using ASV_1, ASV_2, …
#'
#' `r lifecycle::badge("stable")`
#'
#' @inheritParams clean_pq
#'
#' @return A new \code{\link{phyloseq-class}} object with `refseq` slot and new
#'   taxa names
#' @export

add_dna_to_phyloseq <- function(physeq) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(phyloseq::taxa_names(physeq))
  names(dna) <- phyloseq::taxa_names(physeq)
  physeq <- phyloseq::merge_phyloseq(physeq, dna)
  phyloseq::taxa_names(physeq) <-
    paste0("ASV_", seq(phyloseq::ntaxa(physeq)))
  return(physeq)
}
################################################################################



################################################################################
#'  Clean phyloseq object by removing empty samples and taxa
#'
#' @details `r lifecycle::badge("experimental")`
#'
#'  In addition, this function check for discrepancy (and rename) between
#' (i) taxa names in refseq, taxonomy table and otu_table and between
#' (ii) sample names in sam_data and otu_table.
#'
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object obtained
#'   using the `phyloseq` package.
#' @param remove_empty_samples (logical) Do you want to remove samples
#'   without sequences (this is done after removing empty taxa)
#' @param remove_empty_taxa (logical) Do you want to remove taxa
#'   without sequences (this is done before removing empty samples)
#' @param clean_samples_names (logical) Do you want to clean samples names?
#' @param silent (logical) If true, no message are printing.
#' @param verbose (logical) Additional informations in the message
#'   the verbose parameter overwrite the silent parameter.
#' @param force_taxa_as_columns (logical) If true, if the taxa are rows
#'   transpose the otu_table and set taxa_are_rows to false
#' @param force_taxa_as_rows (logical) If true, if the taxa are columns
#'   transpose the otu_table and set taxa_are_rows to true
#' @param reorder_asv (logical) if TRUE the otu_table is ordered by the number of
#'   sequences of ASV (descending order). Default to FALSE.
#'   Need the [microViz](https://github.com/david-barnett/microViz) package.
#' @param rename_asv (logical) if TRUE, ASV are renamed by their position
#'   in the OTU_table (asv_1, asv_2, ...). Default to FALSE. If rename ASV is true,
#'   the ASV names in verbose information can be misleading.
#' @return A new \code{\link{phyloseq-class}} object
#' @export
clean_pq <- function(physeq,
                     remove_empty_samples = TRUE,
                     remove_empty_taxa = TRUE,
                     clean_samples_names = TRUE,
                     silent = FALSE,
                     verbose = FALSE,
                     force_taxa_as_columns = FALSE,
                     force_taxa_as_rows = FALSE,
                     reorder_asv = FALSE,
                     rename_asv = FALSE) {
  verify_pq(physeq)
  if (clean_samples_names) {
    if (!is.null(physeq@refseq)) {
      if (sum(!names(physeq@refseq) %in% taxa_names(physeq)) > 0) {
        names(physeq@refseq) <- taxa_names(physeq)
        if (!silent) {
          message("Change the samples names in refseq slot")
        }
      }
    }
    if (!is.null(physeq@tax_table)) {
      if (sum(!rownames(physeq@tax_table) %in% taxa_names(physeq)) > 0) {
        rownames(physeq@tax_table) <- taxa_names(physeq)
        if (!silent) {
          message("Change the taxa names in tax_table slot")
        }
      }
    }

    if (!is.null(physeq@sam_data)) {
      if (sum(!rownames(physeq@sam_data) %in% sample_names(physeq)) > 0) {
        rownames(physeq@sam_data) <- sample_names(physeq)
        if (!silent) {
          message("Change the samples names in sam_data slot")
        }
      }
    }
  }

  if (reorder_asv) {
    physeq <- microViz::tax_sort(physeq, sum)
  }

  if (rename_asv) {
    taxa_names(physeq) <- paste0("ASV_", seq(1, ntaxa(physeq)))
  }

  if (sum(grepl("^0", "", sample_names(physeq)) > 0) && !silent) {
    message("At least one sample name start with a zero.
    That can be a problem for some phyloseq functions such as
    plot_bar and psmelt.")
  }

  if (force_taxa_as_columns && force_taxa_as_rows) {
    stop("You can't force taxa as column and taxa as row in the same time.")
  }

  if (force_taxa_as_columns && taxa_are_rows(physeq)) {
    otu_table(physeq) <- otu_table(t(as.matrix(unclass(physeq@otu_table))),
      taxa_are_rows = FALSE
    )
    message("Taxa are now in columns.")
  }

  if (force_taxa_as_rows && !taxa_are_rows(physeq)) {
    otu_table(physeq) <- otu_table(t(as.matrix(unclass(physeq@otu_table))),
      taxa_are_rows = TRUE
    )
    message("Taxa are now in rows.")
  }

  new_physeq <- physeq

  if (remove_empty_taxa) {
    if (sum(taxa_sums(new_physeq) == 0) > 0) {
      new_physeq <- subset_taxa(physeq, taxa_sums(physeq) > 0)
    }
  }
  if (remove_empty_samples) {
    if (sum(sample_sums(new_physeq) == 0) > 0) {
      new_physeq <- subset_samples(new_physeq, sample_sums(physeq) > 0)
    }
  }

  if (verbose) {
    message(paste(
      "Cleaning suppress", ntaxa(physeq) - ntaxa(new_physeq), "taxa (",
      paste(taxa_names(physeq)[taxa_sums(physeq) == 0], collapse = " / "), ") and",
      nsamples(physeq) - nsamples(new_physeq),
      "sample(s) (", paste(sample_names(physeq)[sample_sums(physeq) == 0], collapse = " / "), ")."
    ))
  } else if (!silent) {
    message(paste(
      "Cleaning suppress", ntaxa(physeq) - ntaxa(new_physeq), "taxa and",
      nsamples(physeq) - nsamples(new_physeq),
      "samples."
    ))
  }

  return(new_physeq)
}






################################################################################
#' Track the number of reads (= sequences), samples and cluster (e.g. ASV)
#' from various objects including dada-class and derep-class.
#'
#' @details
#' `r lifecycle::badge("maturing")`
#'
#'  * List of fastq and fastg.gz files -> nb of reads and samples
#'  * List of dada-class -> nb of reads, clusters (ASV) and samples
#'  * List of derep-class -> nb of reads, clusters (unique sequences)
#'    and samples
#'  * Matrix of samples x clusters (e.g. `otu_table`) -> nb of reads,
#'    clusters and samples
#'  * Phyloseq-class -> nb of reads, clusters and samples
#'
#' @param list_of_objects (required) a list of objects
#' @param obj_names
#'   A list of names corresponding to the list of objects
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples and empty ASV are discarded
#'   before clustering.
#' @param taxonomy_rank A vector of int. Define the column number of
#'   taxonomic rank `in physeq@tax_table` to compute the number of unique value.
#'   Default is NULL and do not compute values for any taxonomic rank
#' @param ... Other arguments passed on to [clean_pq()] function.
#'
#' @return The number of sequences, clusters (e.g. OTUs, ASVs) and samples for
#'   each object.
#' @export

track_wkflow <- function(
    list_of_objects,
    obj_names = NULL,
    clean_pq = FALSE,
    taxonomy_rank = NULL,
    ...) {
  message("Compute the number of sequences")
  if (!is.null(obj_names)) {
    names(list_of_objects) <- obj_names
  }

  if (clean_pq) {
    for (i in seq_along(list_of_objects)) {
      if (inherits(list_of_objects[[i]], "phyloseq")) {
        list_of_objects[[i]] <- clean_pq(list_of_objects[[i]], ...)
      }
    }
  }

  track_nb_seq_per_obj <-
    pbapply::pblapply(list_of_objects, function(object) {
      message(paste("Start object of class:", class(object), sep = " "))
      if (inherits(object, "phyloseq")) {
        sum(object@otu_table)
      } else if (inherits(object, "matrix")) {
        sum(object, na.rm = TRUE)
      } else if (is.character(object[1]) &&
        length(object[1]) == 1 &&
        file.exists(object[1])) {
        if (summary(file(object[[1]]))$class == "gzfile") {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("zcat ", x, " | grep -c '^+$'", sep = ""),
              intern = TRUE
            ))
          })
        } else if (grepl("\\.fastq$", object[1])) {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("cat ", x, " | grep -c '^+$'", sep = ""),
              intern = TRUE
            ))
          })
        } else {
          stop("Files must be either gzfile or .fastq")
        }
      } else {
        pbapply::pbsapply(object, function(x) {
          sum(dada2::getUniques(x, silence = TRUE))
        })
      }
    })
  track_nb_seq_per_obj <-
    pbapply::pblapply(track_nb_seq_per_obj, sum)

  message("Compute the number of clusters")
  track_nb_cluster_per_obj <-
    pbapply::pblapply(list_of_objects, function(object) {
      message(paste("Start object of class:", class(object), sep = " "))
      if (inherits(object, "phyloseq")) {
        ntaxa(object)
      } else if (inherits(object, "matrix")) {
        ncol(object)
      } else if (inherits(object[[1]], "dada")) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (is.data.frame(object[[1]]) &&
        all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (inherits(object[[1]], "derep")) {
        length(unique(unlist(lapply(object, function(x) {
          names(x$uniques)
        }))))
      } else {
        NA
      }
    })

  message("Compute the number of samples")
  track_nb_sam_per_obj <-
    pbapply::pblapply(list_of_objects, function(object) {
      message(paste("Start object of class:", class(object), sep = " "))
      if (inherits(object, "phyloseq")) {
        nsamples(object)
      } else if (inherits(object, "matrix")) {
        nrow(object)
      } else if (inherits(object[[1]], "dada")) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (is.data.frame(object[[1]]) &&
        all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (inherits(object[[1]], "derep")) {
        length(object)
      } else if (is.character(object[1]) &&
        length(object[1]) == 1 &&
        file.exists(object[1])) {
        length(object)
      } else {
        NA
      }
    })

  if (!is.null(taxonomy_rank)) {
    message("Compute the number of values in taxonomic rank")
    track_nb_tax_value_per_obj <-
      pbapply::pblapply(list_of_objects, function(object) {
        message(paste("Start object of class:", class(object), sep = " "))
        if (inherits(object, "phyloseq")) {
          if (taxa_are_rows(object)) {
            apply(object@tax_table[taxonomy_rank, ], 1, function(x) {
              length(unique(na.omit(x)))
            })
          } else {
            apply(object@tax_table[, taxonomy_rank], 2, function(x) {
              length(unique(na.omit(x)))
            })
          }
        } else {
          rep(NA, length(taxonomy_rank))
        }
      })

    names_taxonomic_rank <-
      pbapply::pblapply(list_of_objects, function(object) {
        message(paste("Start object of class:", class(object), sep = " "))
        if (inherits(object, "phyloseq")) {
          if (taxa_are_rows(object)) {
            rownames(object@tax_table)[taxonomy_rank]
          } else {
            colnames(object@tax_table)[taxonomy_rank]
          }
        }
      })
    track <- plyr::rbind.fill.matrix(
      matrix(ncol = length(list_of_objects), unlist(track_nb_seq_per_obj)),
      matrix(ncol = length(list_of_objects), unlist(track_nb_cluster_per_obj)),
      matrix(ncol = length(list_of_objects), unlist(track_nb_sam_per_obj)),
      matrix(ncol = length(list_of_objects), unlist(track_nb_tax_value_per_obj))
    )

    rownames(track) <- c(
      "nb_sequences",
      "nb_clusters",
      "nb_samples",
      names_taxonomic_rank[[1]]
    )
  } else {
    track <- plyr::rbind.fill.matrix(
      matrix(ncol = length(list_of_objects), unlist(track_nb_seq_per_obj)),
      matrix(ncol = length(list_of_objects), unlist(track_nb_cluster_per_obj)),
      matrix(ncol = length(list_of_objects), unlist(track_nb_sam_per_obj))
    )

    rownames(track) <- c(
      "nb_sequences",
      "nb_clusters",
      "nb_samples"
    )
  }


  track <- as.data.frame(t(track))
  if (!is.null(obj_names)) {
    rownames(track) <- obj_names
  } else {
    rownames(track) <- names(list_of_objects)
  }

  return(track)
}
################################################################################

################################################################################
#' Track the number of reads (= sequences), samples and cluster (e.g. ASV)
#' for each samples
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' More information are available in the manual of the function [track_wkflow()]
#'
#' @param list_pq_obj (required): a list of object passed on to [track_wkflow()]
#' @param ... : other args passed on to [track_wkflow()]
#'
#' @return A list of dataframe. cf [track_wkflow()] for more information
#'
#' @export
#' @md
#'
#' @author Adrien Taudière
#'
#' @examples
#' data(data_fungi)
#' tree_A10_005 <- subset_samples(data_fungi, Tree_name == "A10-005")
#' track_wkflow_samples(tree_A10_005)
track_wkflow_samples <- function(list_pq_obj, ...) {
  if (!inherits(list_pq_obj, "list")) {
    list_pq_obj <- list(list_pq_obj)
  }
  res <- list()
  for (s in sample_names(list_pq_obj[[1]])) {
    list_pq_obj_samples <- lapply(list_pq_obj, select_one_sample, sam_name = s)
    res[[s]] <- track_wkflow(list_pq_obj_samples, ...)
  }
  return(res)
}
################################################################################


################################################################################
#' Recluster sequences of an object of class `physeq`
#' (e.g. OTUs or ASV from dada)
#'
#' `r lifecycle::badge("maturing")`
#'
#' @inheritParams clean_pq
#' @param seq_names You may directly use a character vector of DNA sequences
#'   in place of physeq args. When physeq is set, dna sequences take the value of
#'   `physeq@refseq`
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param method (default: clusterize)
#'   Set the clustering method.
#'   - `clusterize` use the [DECIPHER::Clusterize()] fonction,
#'   - `vsearch` use the vsearch software (https://github.com/torognes/vsearch/)
#'     with arguments `-cluster_fast` and `-strand both`
#' @param vsearchpath path to vsearch
#' @param id (default: 0.97) level of identity to cluster
#' @param tax_adjust See the man page
#'   of [speedyseq::merge_taxa_vec()] for more details.
#' @param vsearch_cluster_method (default: "--cluster_fast") See other possible 
#'   methods in the [vsearch pdf manual](https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch_manual.pdf) (e.g. `--cluster_size` or `--cluster_smallmem`)
#' @param vsearch_args (default : "--strand both") a one length character element defining other parameters to 
#'   passed on to vsearch.
#' @param ... Others arguments path to [DECIPHER::Clusterize()]
#' @details This function use the `speedyseq::merge_taxa_vec` function to
#'   merge taxa into clusters. By default tax_adjust = 1L. See the man page
#'   of [speedyseq::merge_taxa_vec()].
#'
#' @return A new object of class `physeq` or a list of cluster if seq_names
#'   args was used.
#'
#' @references
#'   VSEARCH can be downloaded from
#'   \url{https://github.com/torognes/vsearch}.
#'   More information in the associated publication
#'   \url{https://www.ncbi.nlm.nih.gov/pubmed/27781170}.
#'
#' @export

asv2otu <- function(physeq = NULL,
                    seq_names = NULL,
                    nproc = 1,
                    method = "clusterize",
                    id = 0.97,
                    vsearchpath = "vsearch",
                    tax_adjust = 1,
                    vsearch_cluster_method = "--cluster_fast",
                    vsearch_args = "--strand both",
                    ...) {
  if (inherits(physeq, "phyloseq")) {
    verify_pq(physeq)
    dna <- Biostrings::DNAStringSet(physeq@refseq)
    if (!is.null(seq_names)) {
      stop("You must use either physeq or seq_names args but not both")
    }
  } else if (inherits(seq_names, "character")) {
    dna <- Biostrings::DNAStringSet(seq_names)
  } else {
    stop("You must set the args physeq (object of class phyloseq) or
    seq_names (character vector).")
  }

  if (!method %in% c("clusterize", "vsearch")) {
    stop("Method allows 2 values only : `clusterize` or `vsearch`")
  }

  if (method == "clusterize") {
    ## Find clusters of ASVs to form the new OTUs
    clusters <- DECIPHER::Clusterize(
      dna,
      cutoff = 1 - id,
      # e.g. `cutoff = 0.03` for a 97% OTU
      processors = nproc,
      ...
    )

    if (inherits(physeq, "phyloseq")) {
      new_obj <-
        speedyseq::merge_taxa_vec(physeq,
          clusters$cluster,
          tax_adjust = tax_adjust
        )
    } else if (inherits(seq_names, "character")) {
      new_obj <- clusters
    } else {
      stop("You must set the args physeq (object of class phyloseq) or
    seq_names (character vector).")
    }
  } else if (method == "vsearch") {
    Biostrings::writeXStringSet(dna, "temp.fasta")

    system2(
      vsearchpath,
      paste0(
        paste(" ", vsearch_cluster_method, " temp.fasta ", vsearch_args),
        " -id ",
        id,
        " --centroids cluster.fasta",
        " --uc temp.uc"
      ),
      stdout = TRUE,
      stderr = TRUE
    )

    pack_clusts <- utils::read.table("temp.uc", sep = "\t")
    colnames(pack_clusts) <-
      c(
        "type",
        "cluster",
        "width",
        "identity",
        "strand",
        "6",
        "7",
        "cigarAlignment",
        "query",
        "target"
      )

    clusters <- pack_clusts$cluster
    clusters <-
      tapply(
        clusters, paste(pack_clusts$cluster, pack_clusts$query),
        function(x) {
          x[1]
        }
      )

    if (inherits(physeq, "phyloseq")) {
      new_obj <-
        speedyseq::merge_taxa_vec(physeq,
          clusters,
          tax_adjust = tax_adjust
        )
    } else if (inherits(seq_names, "character")) {
      new_obj <- pack_clusts
    } else {
      stop("You must set the args physeq (object of class phyloseq) or seq_names (character vector).")
    }

    if (file.exists("temp.fasta")) {
      file.remove("temp.fasta")
    }
    if (file.exists("cluster.fasta")) {
      file.remove("cluster.fasta")
    }
    if (file.exists("temp.uc")) {
      file.remove("temp.uc")
    }
  }
  return(new_obj)
}
################################################################################

################################################################################
#' Search for a list of sequence in a fasta file against physeq reference
#'   sequences using [vsearch](https://github.com/torognes/vsearch)
#'
#' `r lifecycle::badge("maturing")`
#'
#' @inheritParams clean_pq
#' @param seq2search (required) path to fasta file
#' @param vsearchpath path to vsearch
#' @param id (default: 0.8) id for the option `--usearch_global` of the vsearch software
#' @param iddef (default: 0) iddef for the option `--usearch_global` of the vsearch software
#' @examples
#' \dontrun{
#' file_dna <- tempfile("dna.fa")
#' seqinr::write.fasta("GCCCATTAGTATTCTAGTGGGCATGCCTGTTCGAGCGTCATTTTCA
#'   ACCCTCAAGCCCCTTATTGCTTGGTGTTGGGAGTTTAGCTGGCTTTATAGCGGTTAACTCCCTAAATATACTGGCG",
#'   file = file_dna, name = "seq1"
#' )
#' res <- vsearch_search_global(data_fungi, file_dna)
#' unlink(file_dna)
#'
#' res[res$identity != "*", ]
#'
#' clean_pq(subset_taxa(data_fungi, res$identity != "*"))
#' }
#' @return A dataframe with uc results (invisible)
#' @export

vsearch_search_global <- function(physeq,
                                  seq2search,
                                  vsearchpath = "vsearch",
                                  id = 0.8,
                                  iddef = 0) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(physeq@refseq)

  Biostrings::writeXStringSet(dna, "temp.fasta")

  system2(
    vsearchpath,
    paste(
      " --usearch_global ",
      here::here("temp.fasta"),
      " --db ",
      here::here(seq2search),
      " --uc",
      " temp.uc",
      " --id ",
      id,
      " --uc_allhits",
      " --strand both",
      " --iddef ",
      iddef,
      sep = ""
    )
  )

  pack_clusts <- utils::read.table("temp.uc", sep = "\t")
  colnames(pack_clusts) <- c(
    "type",
    "cluster",
    "width",
    "identity",
    "strand",
    "6",
    "7",
    "cigarAlignment",
    "query",
    "target"
  )
  unlink("temp.fasta")
  unlink("temp.uc")

  return(invisible(pack_clusts))
}
################################################################################

################################################################################
#' Blast fasta sequences against `refseq` slot of a \code{\link{phyloseq-class}}
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
#'
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
                              list_no_output_query = FALSE) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, "db.fasta")

  system(paste(blastpath,
    "makeblastdb -dbtype nucl -in db.fasta -out dbase",
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
        nrow = length(no_output_query)
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
                     nproc = 1) {
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
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for blast (args -num_threads
#'   for blastn command)
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
                             nproc = 1) {
  blast_tab <- blast_pq(
    physeq = physeq,
    fasta_for_db = fasta_for_db,
    database = database,
    id_cut = id_cut,
    bit_score_cut = bit_score_cut,
    min_cover_cut = min_cover_cut,
    unique_per_seq = TRUE,
    score_filter = TRUE,
    nproc = nproc
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

################################################################################
#' Save phyloseq object in the form of multiple csv tables.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @inheritParams clean_pq
#' @param path a path to the folder to save the phyloseq object
#' @param rdata (logical) does the phyloseq object is also saved in Rdata format?
#' @param one_file (logical) if TRUE, combine all data in one file only
#' @param write_sam_data (logical) does the samples data are add to
#'   the file. Only used if `one_file` is TRUE.
#'   Note that these option result in a lot of NA values.
#' @param sam_data_first (logical) if TRUE, put the sample data at the top of the table
#'   Only used if `one_file` and write_sam_data are both TRUE.
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples are discarded after subsetting ASV
#' @param reorder_asv (logical) if TRUE the otu_table is ordered by the number of
#'   sequences of ASV (descending order). Default to TRUE. Only possible if clean_pq
#'   is set to TRUE. Need the [microViz](https://github.com/david-barnett/microViz) package.
#' @param rename_asv reorder_asv (logical) if TRUE, ASV are renamed by their position
#'   in the OTU_table (asv_1, asv_2, ...). Default to FALSE. Only possible if clean_pq
#'   is set to TRUE.
#' @param quote a logical value (default FALSE) or a numeric vector.
#'   If TRUE, any character or factor columns will be surrounded by
#'   double quotes.  If a numeric vector, its elements are taken
#'   as the indices of columns to quote.  In both cases, row and
#'   column names are quoted if they are written. If FALSE nothing is quoted.
#' @param sep_csv (default tabulation (\t)) separator for column
#' @param ... Other arguments passed on to [utils::write.table()] function.
#' @return One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#'   and if present a phy_tree in Newick format
#' @export
#'
#' @examples
#' \dontrun{
#' write_pq(data_fungi, path = "phyloseq")
#' write_pq(data_fungi, path = "phyloseq", one_file = TRUE)
#' }
#'
write_pq <- function(physeq,
                     path = NULL,
                     rdata = FALSE,
                     one_file = FALSE,
                     write_sam_data = TRUE,
                     sam_data_first = FALSE,
                     clean_pq = TRUE,
                     reorder_asv = FALSE,
                     rename_asv = FALSE,
                     remove_empty_samples = TRUE,
                     remove_empty_taxa = TRUE,
                     clean_samples_names = TRUE,
                     silent = FALSE,
                     verbose = FALSE,
                     quote = FALSE,
                     sep_csv = "\t",
                     ...) {
  verify_pq(physeq)

  physeq <- clean_pq(physeq,
    reorder_asv = reorder_asv,
    rename_asv = rename_asv,
    remove_empty_samples = remove_empty_samples,
    remove_empty_taxa = remove_empty_taxa,
    clean_samples_names = clean_samples_names,
    silent = silent,
    verbose = verbose
  )

  if (!dir.exists(path)) {
    dir.create(file.path(path), recursive = TRUE)
  }
  if (one_file) {
    if (!is.null(physeq@refseq) && !is.null(physeq@otu_table) && !is.null(physeq@tax_table)) {
      if (!taxa_are_rows(physeq)) {
        otu_table(physeq) <- otu_table(t(as.matrix(unclass(physeq@otu_table))),
          taxa_are_rows = TRUE
        )
      }
      df_physeq_interm <- cbind(
        physeq@otu_table,
        physeq@tax_table,
        as.vector(physeq@refseq)
      )
      colnames(df_physeq_interm) <- c(sample_names(physeq), colnames(physeq@tax_table), "Reference Sequences")

      df_physeq_interm <- as.data.frame(df_physeq_interm)

      if (write_sam_data) {
        sam_data <- data.frame(t(data.frame(unclass(physeq@sam_data))))
        colnames(sam_data) <- sample_names(physeq)
        if (sam_data_first) {
          df_physeq <- dplyr::full_join(sam_data, df_physeq_interm)
          rownames(df_physeq) <- c(rownames(sam_data), rownames(df_physeq_interm))
        } else {
          df_physeq <- dplyr::full_join(df_physeq_interm, sam_data)
          rownames(df_physeq) <- c(rownames(df_physeq_interm), rownames(sam_data))
        }
      } else {
        df_physeq <- df_physeq_interm
      }
      utils::write.table(
        df_physeq,
        paste0(path, "/ASV_table_allInOne.csv"),
        quote = quote,
        ...
      )
    } else if (!is.null(physeq@otu_table) && !is.null(physeq@tax_table)) {
      if (!taxa_are_rows(physeq)) {
        otu_table(physeq) <- otu_table(t(as.matrix(unclass(physeq@otu_table))),
          taxa_are_rows = TRUE
        )
      }
      df_physeq_interm <- cbind(
        physeq@otu_table,
        physeq@tax_table
      )
      colnames(df_physeq_interm) <- c(sample_names(physeq), colnames(physeq@tax_table), "Reference Sequences")

      df_physeq_interm <- as.data.frame(df_physeq_interm)

      if (write_sam_data) {
        sam_data <- data.frame(t(data.frame(unclass(physeq@sam_data))))
        colnames(sam_data) <- sample_names(physeq)
        if (sam_data_first) {
          df_physeq <- dplyr::full_join(sam_data, df_physeq_interm)
          rownames(df_physeq) <- c(rownames(sam_data), rownames(df_physeq_interm))
        } else {
          df_physeq <- dplyr::full_join(df_physeq_interm, sam_data)
          rownames(df_physeq) <- c(rownames(df_physeq_interm), rownames(sam_data))
        }
      }
      utils::write.table(
        df_physeq,
        paste0(path, "/ASV_table_allInOne.csv"),
        quote = quote,
        sep = sep_csv,
        ...
      )
    }
  } else {
    if (!is.null(physeq@otu_table)) {
      utils::write.table(
        physeq@otu_table, 
        paste0(path, "/otu_table.csv"),
        quote = quote,
        sep = sep_csv,
        ...
      )
    }
    if (!is.null(physeq@refseq)) {
      utils::write.table(
        physeq@refseq,
        paste0(path, "/refseq.csv"),
        quote = quote,
        sep = sep_csv,
        ...
      )
    }
    if (!is.null(physeq@tax_table)) {
      utils::write.table(
        physeq@tax_table, 
        paste0(path, "/tax_table.csv"),
        quote = quote,
        sep = sep_csv,
        col.names = NA,
        ...
      )
    }
    if (!is.null(physeq@sam_data)) {
      utils::write.table(
        as.matrix(physeq@sam_data), 
        paste0(path, "/sam_data.csv"),
        quote = quote,
        sep = sep_csv,
        col.names = NA,
        ...
      )
    }
  }
  if (!is.null(physeq@phy_tree)) {
    ape::write.tree(physeq@phy_tree, paste(path, "/phy_tree.txt", sep = ""))
  }
  if (rdata) {
    save(physeq, file = paste(path, "/physeq.RData", sep = ""))
  }
}
################################################################################

################################################################################
#' Read phyloseq object from multiple csv tables and a phylogenetic tree
#' in Newick format.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @param path (required) a path to the folder to read the phyloseq object
#' @param taxa_are_rows (required, default to FALSE) see ?phyloseq for details
#' @param sam_names The name of the variable (column) in sam_data.csv to rename
#'   samples. Note that if you use [write_phyloseq()] function to save your
#'   physeq object, you may use sam_names = "X" to rename the samples names
#'   as before.
#' @param sep_csv (default tabulation (\t)) separator for column
#' @param ... Other arguments passed on to [utils::write.table()] function.
#' @return One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#' and if present a phy_tree in Newick format. At least the otu_table.csv need to be present.
#' @export
#'
#' @examples
#' \dontrun{
#' read_pq(path = "phyloseq_data")
#' }
#'
read_pq <- function(path = NULL, taxa_are_rows = FALSE, sam_names = NULL, sep_csv = "\t", ...) {
  
if (file.exists(paste0(path, "/otu_table.csv"))) {
   if (taxa_are_rows) {
    otu_table_csv <- as.matrix(utils::read.table(paste0(path, "/otu_table.csv"), sep = sep_csv))
    samp_names <- colnames(otu_table_csv)
    otu_table_csv <- apply(otu_table_csv, 2, as.numeric)
    table_otu <- otu_table(otu_table_csv, taxa_are_rows = TRUE)
    sample_names(table_otu) <- samp_names
    physeq <- phyloseq(table_otu)
  } else {
    otu_table_csv <- as.matrix(utils::read.table(paste0(path, "/otu_table.csv"), sep = sep_csv))
    samp_names <- rownames(otu_table_csv)
    otu_table_csv <- apply(otu_table_csv, 2, as.numeric)
    rownames(otu_table_csv) <-  samp_names
    physeq <- phyloseq(otu_table(otu_table_csv, taxa_are_rows = FALSE))
  }
}
if (file.exists(paste0(path, "/refseq.csv"))) {
  dna <- Biostrings::DNAStringSet(utils::read.table(paste0(path, "/refseq.csv"),  sep = sep_csv, row.names = NULL)[, 2])
  names(dna) <- utils::read.table(paste0(path, "/refseq.csv"), sep = sep_csv, row.names = NULL)[, 1]
  physeq <- phyloseq::merge_phyloseq(physeq, refseq(dna))
}
if (file.exists(paste0(path, "/tax_table.csv"))) {
  tax_table_csv <- utils::read.table(paste0(path, "/tax_table.csv"), sep = sep_csv)
  rownames(tax_table_csv) <- tax_table_csv[, 1]
  tax_table_csv <- as.matrix(tax_table_csv[, -1])
  physeq <- phyloseq::merge_phyloseq(physeq, tax_table(tax_table_csv))
}
if (file.exists(paste0(path, "/sam_data.csv"))) {
  sam_data_csv <- utils::read.table(paste0(path, "/sam_data.csv"), sep = sep_csv)
  rownames(sam_data_csv) <- sam_data_csv[,1]
  physeq <- phyloseq::merge_phyloseq(physeq, sample_data(sam_data_csv))
}

if (!is.null(physeq@phy_tree)) {
  tree <- ape::read.tree(paste0(path, "/phy_tree.txt"))
  physeq <- phyloseq::merge_phyloseq(physeq, phy_tree(tree))
}
if (!is.null(sam_names)) {
  sample_names(physeq) <- unclass(physeq@sam_data[, sam_names])[[1]]
}


  return(physeq)
}

################################################################################

################################################################################
#' Lulu reclustering of class `physeq`
#'
#' @details
#' `r lifecycle::badge("experimental")`
#'
#' See https://www.nature.com/articles/s41467-017-01312-x for more information
#'  on the method.

#'
#' @inheritParams clean_pq
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param id (default: 0.84) id for --usearch_global.
#' @param vsearchpath path to vsearch.
#' @param verbose (logical) if true, print some additional messages.
#' @param clean_pq (logical) if true, empty samples and empty ASV are discarded
#'   before clustering.
#'
#' @return a list of for object
#' - "new_physeq": The new phyloseq object (class physeq)
#' - "discrepancy_vector": A vector of discrepancy showing for each taxonomic
#'   level the proportion of identic value before and after lulu reclustering.
#'   A value of 0.6 stands for 60% of ASV before re-clusering have
#'   identical value after re-clustering. In other word, 40% of ASV are assigned
#'   to a different taxonomic
#'   value. NA value are not counted as discrepancy.
#' - "res_lulu": A list of the result from the lulu function
#' - "merged_ASV": the data.frame used to merged ASV
#'
#' @importFrom stats na.exclude
#' @export
#' @examples
#' \dontrun{
#' data(data_fungi_sp_known)
#' lulu_pq(data_fungi_sp_known)
#' }
#' @author Tobias Guldberg Frøslev \email{tobiasgf@snm.ku.dk}
#'   & Adrien Taudière \email{adrien.taudiere@@zaclys.net}
#' @details
#' The version of LULU is a fork of Adrien Taudière (\url{https://github.com/adrientaudiere/lulu})
#'  from \url{https://github.com/tobiasgf/lulu}
#' @references
#' - LULU : \url{https://github.com/adrientaudiere/lulu}
#'  forked from \url{https://github.com/tobiasgf/lulu}.
#' - VSEARCH can be downloaded from
#'  \url{https://github.com/torognes/vsearch}.

lulu_pq <- function(physeq,
                    nproc = 1,
                    id = 0.84,
                    vsearchpath = "vsearch",
                    verbose = FALSE,
                    clean_pq = FALSE) {
  verify_pq(physeq)
  if (clean_pq) {
    physeq <- clean_pq(physeq)
  }

  message("Start Vsearch usearch_global")
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, "temp.fasta")
  system2(
    vsearchpath,
    paste(
      " --usearch_global temp.fasta --db temp.fasta --self --iddef 1",
      " -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10",
      " -id ",
      id,
      "  --userout match_list.txt",
      sep = ""
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  match_list <- utils::read.csv(file = "match_list.txt", sep = "\t")

  if (!requireNamespace("lulu")) {
    requireNamespace(devtools)
    install_github("adrientaudiere/lulu")
  }

  message("Lulu algorithm")
  res_lulu <-
    lulu::lulu(data.frame(t(physeq@otu_table)), match_list)

  if (file.exists("temp.fasta")) {
    file.remove("temp.fasta")
  }
  if (file.exists("cluster.fasta")) {
    file.remove("cluster.fasta")
  }
  if (file.exists("temp.uc")) {
    file.remove("temp.uc")
  }
  if (file.exists("match_list.txt")) {
    file.remove("match_list.txt")
  }

  merged <- res_lulu$otu_map[res_lulu$otu_map$curated == "merged", ]
  merged <- merged[rownames(merged) != merged$parent_id, ]

  test_vector <- c()
  for (tax_rank in colnames(physeq@tax_table)) {
    test <- physeq@tax_table[rownames(merged), tax_rank] == physeq@tax_table[merged$parent_id, tax_rank]
    test_vector <- c(test_vector, sum(test, na.rm = TRUE) / length(stats::na.exclude(test)))
  }

  names(test_vector) <- colnames(physeq@tax_table)

  new_physeq <- prune_taxa(taxa_names(physeq) %in% rownames(res_lulu$curated_table), physeq)
  new_physeq@otu_table <- otu_table(t(res_lulu$curated_table), taxa_are_rows = FALSE)
  sample_names(new_physeq) <- sample_names(physeq)

  if (verbose) {
    message(paste("The number of taxa decrease from ", ntaxa(physeq), " to ", ntaxa(new_physeq), ".", sep = ""))
    message("See the discrepancy_vector to verify the degree of discrepancy in taxonomy due to lulu re-clustering.")
  }
  return(list(
    "new_physeq" = new_physeq,
    "discrepancy_vector" = test_vector,
    "res_lulu" = res_lulu,
    "merged_ASV" = merged
  ))
}
################################################################################

################################################################################
#' Verify the validity of a phyloseq object
#'
#' @details
#' `r lifecycle::badge("maturing")`
#'
#' Mostly for internal use in MiscMetabar functions.
#'
#' @inheritParams clean_pq
#'
#' @return Nothing if the phyloseq object is valid. An error in the other case.
#' @export
#'
verify_pq <- function(physeq) {
  if (!methods::validObject(physeq) || !inherits(physeq, "phyloseq")) {
    stop("The physeq argument is not a valid phyloseq object.")
  }
}
################################################################################


################################################################################
#' Subset samples using a conditional boolean vector.
#'
#' @details
#' `r lifecycle::badge("experimental")`
#'
#' The main objective of this function is to complete the [phyloseq::subset_samples()]
#' function by propose a more easy (but more prone to error) way of subset_samples.
#' It replace the subsetting expression which used the name of the variable
#' in the sam_data by a boolean vector.
#'
#' Warnings: you must verify the result of this function as the
#' boolean condition must match the order of samples in the `sam_data`
#' slot.
#'
#' @inheritParams clean_pq
#' @param condition A boolean vector to subset samples. Length must fit
#'   the number of samples
#'
#' @return a new phyloseq object
#' @export
#'
subset_samples_pq <- function(physeq, condition) {
  if (length(condition) != nsamples(physeq)) {
    stop("Length of condition is different from the number of samples.")
  }
  if (is.null(sample_data(physeq))) {
    cat("Nothing subset. No sample_data in physeq.\n")
    return(physeq)
  } else {
    oldDF <- as(sample_data(physeq), "data.frame")
    newDF <- oldDF[condition, ]
    if (class(physeq) == "sample_data") {
      return(sample_data(newDF))
    } else {
      sample_data(physeq) <- sample_data(newDF)
      return(physeq)
    }
  }
}
################################################################################

################################################################################
#' Subset taxa using a conditional named boolean vector.
#'
#' @details
#' `r lifecycle::badge("experimental")`
#'
#' The main objective of this function is to complete the [phyloseq::subset_taxa()]
#' function by propose a more easy way of subset_taxa using a named boolean vector.
#' Names must match taxa_names.
#'
#'
#' @inheritParams clean_pq
#' @param condition A named boolean vector to subset taxa. Length must fit
#'   the number of taxa and names must match taxa_names
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples are discarded after subsetting ASV
#' @param verbose (logical) Informations are printed

#' @return a new phyloseq object
#' @export
#'
subset_taxa_pq <- function(physeq, condition, verbose = TRUE, clean_pq = TRUE) {
  if (!sum(names(condition) %in% taxa_names(physeq)) == length(condition)) {
    stop(paste(
      "Some names in condition do not fit taxa_names of physeq : ",
      paste(names(condition)[!names(condition) %in% taxa_names(physeq)], collapse = "/")
    ))
  }

  new_physeq <- physeq

  if (!taxa_are_rows(new_physeq)) {
    new_physeq@otu_table <- otu_table(t(new_physeq@otu_table), taxa_are_rows = TRUE)
    taxa_are_rows(new_physeq) <- TRUE
  }

  cond <- condition[match(taxa_names(new_physeq), names(condition))]
  cond[is.na(cond)] <- FALSE

  oldMA <- as(otu_table(new_physeq), "matrix")
  newMA <- oldMA[cond, ]

  if (!is.matrix(newMA)) {
    newMA <- as.matrix(newMA)
    new_otu_table <- otu_table(newMA, taxa_are_rows = TRUE)
    sample_names(new_otu_table) <- sample_names(new_physeq)
  } else {
    new_otu_table <- otu_table(newMA, taxa_are_rows = TRUE)
  }

  otu_table(new_physeq) <- new_otu_table

  if (clean_pq) {
    new_physeq <- clean_pq(new_physeq, verbose = TRUE)
  }

  if (verbose) {
    message(paste("Number of non-matching ASV", sum(is.na(match(taxa_names(physeq), names(condition))))))
    message(paste("Number of matching ASV", sum(!is.na(match(taxa_names(physeq), names(condition))))))
    message(paste("Number of filtered-out ASV", ntaxa(physeq) - ntaxa(new_physeq)))
    message(paste("Number of kept ASV", ntaxa(new_physeq)))
    message(paste("Number of kept samples", nsamples(new_physeq)))
  }

  return(new_physeq)
}
################################################################################


################################################################################
#' Select one sample from a physeq object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Mostly for internal used, for example in function [track_wkflow_samples()].
#'
#' @inheritParams clean_pq
#' @param sam_name (required) The sample name to select
#' @param silent (logical) If true, no message are printing.
#' @return a physeq object with one sample
#'
#' @export
#' @md
#'
#' @author Adrien Taudière
#'
#' @examples
#' data(data_fungi)
#' A8_005 <- select_one_sample(data_fungi, "A8-005_S4_MERGED.fastq.gz")
#' A8_005
select_one_sample <- function(physeq, sam_name, silent = FALSE) {
  if (sum(sample_names(physeq) %in% sam_name) == 0) {
    stop(paste0("The sample ", sam_name, " is not present in the names of samples of your phyloseq physeq object. You may use the sample_names() function."))
  }
  cl_sam <- clean_pq(subset_samples_pq(physeq, sample_names(physeq) == sam_name), silent = TRUE)

  if (!silent) {
    message(paste0("You select 1 of ", nsamples(physeq), " samples and conserved ", ntaxa(cl_sam), " out of ", ntaxa(physeq), " taxa represented by ", sum(cl_sam@otu_table), " sequences (out of ", sum(physeq@otu_table), " sequences [", perc(sum(cl_sam@otu_table), sum(physeq@otu_table)), "%])"))
  }

  return(cl_sam)
}
################################################################################
