if (getRversion() >= "2.15.1") {
  utils::globalVariables(".")
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
#' @description `r lifecycle::badge("experimental")`
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
#' @param rename_asv (logical) if TRUE, ASV are renamed by their position
#'   in the OTU_table (asv_1, asv_2, ...). Default to FALSE. If rename ASV is true,
#'   the ASV names in verbose information can be misleading.
#' @param simplify_taxo (logical) if TRUE, correct the taxonomy_table using the
#'   `MiscMetabar::simplify_taxo()` function
#' @return A new \code{\link{phyloseq-class}} object
#' @export
clean_pq <- function(
        physeq,
        remove_empty_samples = TRUE,
        remove_empty_taxa = TRUE,
        clean_samples_names = TRUE,
        silent = FALSE,
        verbose = FALSE,
        force_taxa_as_columns = FALSE,
        force_taxa_as_rows = FALSE,
        reorder_asv = FALSE,
        rename_asv = FALSE,
        simplify_taxo = FALSE) {
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

  verify_pq(physeq)

  if (reorder_asv) {
    physeq <- reorder_taxa_pq(
      physeq,
      taxa_names(physeq)[order(taxa_sums(physeq), decreasing = TRUE)]
    )
  }

  if (rename_asv) {
    taxa_names(physeq) <- paste0("ASV_", seq(1, ntaxa(physeq)))
  }

  if (sum(grepl("^0", sample_names(physeq)) > 0) && !silent) {
    message(
      "At least one sample name start with a zero.
    That can be a problem for some phyloseq functions such as
    plot_bar and psmelt."
    )
  }

  if (force_taxa_as_columns && force_taxa_as_rows) {
    stop("You can't force taxa as column and taxa as row in the same time.")
  }

  if (force_taxa_as_columns && taxa_are_rows(physeq)) {
    otu_table(physeq) <-
      otu_table(
        t(as.matrix(unclass(
          physeq@otu_table
        ))),
        taxa_are_rows = FALSE
      )
    message("Taxa are now in columns.")
  }

  if (force_taxa_as_rows && !taxa_are_rows(physeq)) {
    otu_table(physeq) <-
      otu_table(
        t(as.matrix(unclass(
          physeq@otu_table
        ))),
        taxa_are_rows = TRUE
      )
    message("Taxa are now in rows.")
  }

  if (simplify_taxo) {
    physeq <- simplify_taxo(physeq)
  }

  new_physeq <- physeq

  if (remove_empty_taxa) {
    if (sum(taxa_sums(new_physeq) != 0) > 0) {
      new_physeq <- subset_taxa(physeq, taxa_sums(physeq) > 0)
    }
  }
  if (remove_empty_samples) {
    if (sum(sample_sums(new_physeq) != 0) > 0) {
      new_physeq <- subset_samples(new_physeq, sample_sums(physeq) > 0)
    }
  }

  if (verbose) {
    message(
      paste(
        "Cleaning suppress",
        ntaxa(physeq) - ntaxa(new_physeq),
        "taxa (",
        paste(taxa_names(physeq)[taxa_sums(physeq) == 0], collapse = " / "),
        ") and",
        nsamples(physeq) - nsamples(new_physeq),
        "sample(s) (",
        paste(sample_names(physeq)[sample_sums(physeq) == 0], collapse = " / "),
        ")."
      )
    )
  } else if (!silent) {
    message(
      paste(
        "Cleaning suppress",
        ntaxa(physeq) - ntaxa(new_physeq),
        "taxa and",
        nsamples(physeq) - nsamples(new_physeq),
        "samples."
      )
    )
  }

  verify_pq(new_physeq)
  return(new_physeq)
}






################################################################################
#' Track the number of reads (= sequences), samples and cluster (e.g. ASV)
#' from various objects including dada-class and derep-class.
#'
#' @description
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
      } else if (inherits(object, "derep")) {
        sum(object$uniques)
      } else if (inherits(object, "dada")) {
        sum(dada2::getUniques(object))
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
      } else if (inherits(object, "dada")) {
        length(object$sequence)
      } else if (inherits(object[[1]], "dada")) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (is.data.frame(object[[1]]) &&
        all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (inherits(object, "derep")) {
        length(unique(names(object$uniques)))
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
      } else if (inherits(object, "dada")) {
        1
      } else if (inherits(object[[1]], "dada")) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (is.data.frame(object[[1]]) &&
        all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (inherits(object, "derep")) {
        1
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
              length(unique(stats::na.omit(x)))
            })
          } else {
            apply(object@tax_table[, taxonomy_rank], 2, function(x) {
              length(unique(stats::na.omit(x)))
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
      matrix(
        ncol = length(list_of_objects),
        unlist(track_nb_cluster_per_obj)
      ),
      matrix(ncol = length(list_of_objects), unlist(track_nb_sam_per_obj)),
      matrix(
        ncol = length(list_of_objects),
        unlist(track_nb_tax_value_per_obj)
      )
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
      matrix(
        ncol = length(list_of_objects),
        unlist(track_nb_cluster_per_obj)
      ),
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
#' for each samples.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Contrary to [track_wkflow()], only phyloseq object are possible.
#' More information are available in the manual of the function [track_wkflow()]
#'
#' @param list_pq_obj (required): a list of object passed on to [track_wkflow()]
#'   Only phyloseq object will return value because information of sample is needed
#' @param ... : other args passed on to [track_wkflow()]
#'
#' @return A list of dataframe. cf [track_wkflow()] for more information
#'
#' @export
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
  if (sum(!unlist(lapply(list_pq_obj, inherits, "phyloseq"))) != 0) {
    stop("At least one object in your list_pq_obj is not a phyloseq obj.")
  }
  res <- list()
  sam_names <- unique(unlist(lapply(list_pq_obj, sample_names)))
  for (s in sam_names) {
    list_pq_obj_samples <-
      lapply(list_pq_obj, function(physeq) {
        if (sum(sample_names(physeq) %in% s) == 1) {
          select_one_sample(physeq, sam_name = s)
        } else {
          matrix(0, nrow = 0, ncol = 0)
        }
      })
    res[[s]] <- track_wkflow(list_pq_obj_samples) # ,...)
  }
  return(res)
}
###########################################################################


################################################################################
#' Recluster sequences of an object of class `physeq`
#' (e.g. OTUs or ASV from dada)
#'
#' @description
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
#'     with arguments `--cluster_size` by default (see args `vsearch_cluster_method`)
#'     and `-strand both` (see args `vsearch_args`)
#' @param vsearchpath path to vsearch
#' @param id (default: 0.97) level of identity to cluster
#' @param tax_adjust See the man page
#'   of [merge_taxa_vec()] for more details.
#'   To conserved the taxonomic rank of the most abundant ASV,
#'   set tax_adjust to 0
#' @param vsearch_cluster_method (default: "--cluster_size) See other possible
#'   methods in the [vsearch pdf manual](https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch_manual.pdf) (e.g. `--cluster_size` or `--cluster_smallmem`)
#'   - `--cluster_fast` : Clusterize the fasta sequences in filename, automatically sort by decreasing sequence length beforehand.
#'   - `--cluster_size` : Clusterize the fasta sequences in filename, automatically sort by decreasing sequence abundance beforehand.
#'   - `--cluster_smallmem` : Clusterize the fasta sequences in filename without automatically modifying their order beforehand. Sequence are expected to be sorted by decreasing sequence length, unless *--usersort* is used
#' @param vsearch_args (default : "--strand both") a one length character element defining other parameters to
#'   passed on to vsearch.
#' @param  keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - temp.fasta (refseq in fasta)
#'   - cluster.fasta (centroid if method = "vsearch")
#'   - temp.uc (clusters if method = "vsearch")
#' @param ... Others arguments path to [DECIPHER::Clusterize()]
#' @details This function use the `merge_taxa_vec` function to
#'   merge taxa into clusters. By default tax_adjust = 1L. See the man page
#'   of [merge_taxa_vec()].
#'
#' @return A new object of class `physeq` or a list of cluster if seq_names
#'   args was used.
#'
#' @references
#'   VSEARCH can be downloaded from
#'   \url{https://github.com/torognes/vsearch}.
#'   More information in the associated publication
#'   \url{https://pubmed.ncbi.nlm.nih.gov/27781170}.
#'
#' @export

asv2otu <- function(
        physeq = NULL,
        seq_names = NULL,
        nproc = 1,
        method = "clusterize",
        id = 0.97,
        vsearchpath = "vsearch",
        tax_adjust = 0,
        vsearch_cluster_method = "--cluster_size",
        vsearch_args = "--strand both",
        keep_temporary_files = FALSE,
        ...) {
  if (inherits(physeq, "phyloseq")) {
    verify_pq(physeq)
    if (is.null(physeq@refseq)) {
      stop("The phyloseq object do not contain a @refseq slot")
    }
    dna <- Biostrings::DNAStringSet(physeq@refseq)
    if (!is.null(seq_names)) {
      stop("You must use either physeq or seq_names args but not both")
    }
  } else if (inherits(seq_names, "character")) {
    dna <- Biostrings::DNAStringSet(seq_names)
  } else {
    stop(
      "You must set the args physeq (object of class phyloseq) or
    seq_names (character vector)."
    )
  }

  if (!method %in% c("clusterize", "vsearch")) {
    stop("Method allows 2 values only : `clusterize` or `vsearch`")
  }

  if (method == "clusterize") {
    ## Find clusters of ASVs to form the new OTUs
    clusters <- DECIPHER::Clusterize(dna,
      cutoff = 1 - id,
      # e.g. `cutoff = 0.03` for a 97% OTU
      processors = nproc,
      ...
    )

    if (inherits(physeq, "phyloseq")) {
      new_obj <-
        merge_taxa_vec(physeq,
          clusters$cluster,
          tax_adjust = tax_adjust
        )
    } else if (inherits(seq_names, "character")) {
      new_obj <- clusters
    } else {
      stop(
        "You must set the args physeq (object of class phyloseq) or
    seq_names (character vector)."
      )
    }
  } else if (method == "vsearch") {
    Biostrings::writeXStringSet(dna, paste0(tempdir(), "/", "temp.fasta"))

    system2(
      vsearchpath,
      paste0(
        paste0(
          " ",
          vsearch_cluster_method,
          " ",
          paste0(tempdir(), "/", "temp.fasta"),
          " ",
          vsearch_args
        ),
        " -id ",
        id,
        " --centroids ",
        paste0(tempdir(), "/", "cluster.fasta"),
        " --uc ",
        paste0(tempdir(), "/", "temp.uc")
      ),
      stdout = TRUE,
      stderr = TRUE
    )

    pack_clusts <-
      utils::read.table(paste0(tempdir(), "/", "temp.uc"), sep = "\t")
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

    clusters <- pack_clusts$cluster[pack_clusts$type != "C"]
    names(clusters) <- pack_clusts$query[pack_clusts$type != "C"]
    clusters <- clusters[match(taxa_names(physeq), names(clusters))]

    if (inherits(physeq, "phyloseq")) {
      new_obj <-
        merge_taxa_vec(physeq,
          clusters,
          tax_adjust = tax_adjust
        )
    } else if (inherits(seq_names, "character")) {
      new_obj <- pack_clusts
    } else {
      stop(
        "You must set the args physeq (object of class phyloseq) or seq_names (character vector)."
      )
    }

    if (file.exists(paste0(tempdir(), "/", "temp.fasta")) &&
      !keep_temporary_files) {
      unlink(paste0(tempdir(), "/", "temp.fasta"))
    }
    if (file.exists(paste0(tempdir(), "/", "cluster.fasta")) &&
      !keep_temporary_files) {
      unlink(paste0(tempdir(), "/", "cluster.fasta"))
    }
    if (file.exists(paste0(tempdir(), "/", "temp.uc")) &&
      !keep_temporary_files) {
      unlink(paste0(tempdir(), "/", "temp.uc"))
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
#' @param seq2search (required if path_to_fasta is NULL) Either (i) a DNAstringSet object
#'   or (ii) a character vector that will be convert to DNAstringSet using
#'   [Biostrings::DNAStringSet()]
#' @param path_to_fasta (required if seq2search is NULL) a path to fasta file if seq2search is est to NULL.
#' @param vsearchpath (default: vsearch) path to vsearch
#' @param id (default: 0.8) id for the option `--usearch_global` of the vsearch software
#' @param iddef (default: 0) iddef for the option `--usearch_global` of the vsearch software
#' @param  keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - temp.fasta (refseq in fasta)
#'   - cluster.fasta (centroid)
#'   - temp.uc (clusters)
#' @examples
#' \dontrun{
#' file_dna <- tempfile("dna.fa")
#' seqinr::write.fasta("GCCCATTAGTATTCTAGTGGGCATGCCTGTTCGAGCGTCATTTTCA
#'   ACCCTCAAGCCCCTTATTGCTTGGTGTTGGGAGTTTAGCTGGCTTTATAGCGGTTAACTCCCTAAATATACTGGCG",
#'   file = file_dna, names = "seq1"
#' )
#' data(data_fungi)
#' res <- vs_search_global(data_fungi, file_dna)
#' unlink(file_dna)
#'
#' res[res$identity != "*", ]
#'
#' clean_pq(subset_taxa(data_fungi, res$identity != "*"))
#' }
#' @return A dataframe with uc results (invisible)
#' @export
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make [vsearch](https://github.com/torognes/vsearch).

vs_search_global <- function(
        physeq,
        seq2search = NULL,
        path_to_fasta = NULL,
        vsearchpath = "vsearch",
        id = 0.8,
        iddef = 0,
        keep_temporary_files = FALSE) {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, paste0(tempdir(), "/", "temp.fasta"))

  if (is.null(seq2search) && is.null(path_to_fasta)) {
    stop("You must fill either seq2search or path_to_fasta argument.")
  }

  if (!is.null(seq2search) && !is.null(path_to_fasta)) {
    stop("You must set either seq2search or path_to_fasta but not both.")
  }
  if (!is.null(seq2search)) {
    if (inherits(seq2search, "character")) {
      seq2search <- Biostrings::DNAStringSet(seq2search)
    }
    Biostrings::writeXStringSet(seq2search, paste0(tempdir(), "seq2search.fasta"))
    seq2search <- paste0(tempdir(), "seq2search.fasta")
  } else if (!is.null(path_to_fasta)) {
    dna <- Biostrings::readDNAStringSet(path_to_fasta)
    Biostrings::writeXStringSet(dna, paste0(tempdir(), "seq2search.fasta"))
    seq2search <- paste0(tempdir(), "seq2search.fasta")
  }

  system2(
    vsearchpath,
    paste(
      " --usearch_global ",
      paste0(tempdir(), "/", "temp.fasta"),
      " --db ",
      seq2search,
      " --uc ",
      paste0(tempdir(), "/", "temp.uc"),
      " --id ",
      id,
      " --uc_allhits",
      " --strand both",
      " --iddef ",
      iddef,
      sep = ""
    )
  )

  pack_clusts <-
    utils::read.table(paste0(tempdir(), "/", "temp.uc"), sep = "\t")
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

  if (!keep_temporary_files) {
    if(file.exists(paste0(tempdir(), "temp.fasta"))){
      unlink(paste0(tempdir(), "temp.fasta"))
    }
    if(file.exists(paste0(tempdir(), "temp.uc"))){
      unlink(paste0(tempdir(), "temp.uc"))
    }
    if(file.exists(paste0(tempdir(), "seq2search.fasta"))){
      unlink(paste0(tempdir(), "seq2search.fasta"))
    }
  } else {
    message(paste0("Temporary files are located at ", tempdir()))
  }

  return(invisible(pack_clusts))
}
################################################################################



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
#'   is set to TRUE.
#' @param rename_asv reorder_asv (logical) if TRUE, ASV are renamed by their position
#'   in the OTU_table (asv_1, asv_2, ...). Default to FALSE. Only possible if clean_pq
#'   is set to TRUE.
#' @param quote a logical value (default FALSE) or a numeric vector.
#'   If TRUE, any character or factor columns will be surrounded by
#'   double quotes.  If a numeric vector, its elements are taken
#'   as the indices of columns to quote.  In both cases, row and
#'   column names are quoted if they are written. If FALSE nothing is quoted.
#' @param sep_csv (default tabulation) separator for column
#' @param ... Other arguments passed on to [utils::write.table()] function.
#' @return Build a folder (path) containing one to four csv tables
#'   (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#'   and if present a phy_tree in Newick format
#' @export
#' @author Adrien Taudière
#' @examples
#' data(data_fungi)
#' # write_pq(data_fungi, path = "phyloseq")
#' # write_pq(data_fungi, path = "phyloseq", one_file = TRUE)
#' @seealso [MiscMetabar::save_pq()]

write_pq <- function(
        physeq,
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

  physeq <- clean_pq(
    physeq,
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
    if (!is.null(physeq@refseq) &&
      !is.null(physeq@otu_table) && !is.null(physeq@tax_table)) {
      if (!taxa_are_rows(physeq)) {
        otu_table(physeq) <-
          otu_table(
            t(as.matrix(unclass(
              physeq@otu_table
            ))),
            taxa_are_rows = TRUE
          )
      }
      df_physeq_interm <- cbind(
        physeq@otu_table,
        physeq@tax_table,
        as.vector(physeq@refseq)
      )
      colnames(df_physeq_interm) <-
        c(
          sample_names(physeq),
          colnames(physeq@tax_table),
          "Reference Sequences"
        )

      df_physeq_interm <- as.data.frame(df_physeq_interm)

      if (write_sam_data) {
        sam_data <- data.frame(t(data.frame(unclass(
          physeq@sam_data
        ))))
        colnames(sam_data) <- sample_names(physeq)
        if (sam_data_first) {
          df_physeq <- dplyr::full_join(sam_data, df_physeq_interm)
          rownames(df_physeq) <-
            c(rownames(sam_data), rownames(df_physeq_interm))
        } else {
          df_physeq <- dplyr::full_join(df_physeq_interm, sam_data)
          rownames(df_physeq) <-
            c(rownames(df_physeq_interm), rownames(sam_data))
        }
      } else {
        df_physeq <- df_physeq_interm
      }
      utils::write.table(
        df_physeq,
        paste0(path, "/ASV_table_allInOne.csv"),
        quote = quote,
        sep = sep_csv,
        ...
      )
    } else if (!is.null(physeq@otu_table) &&
      !is.null(physeq@tax_table)) {
      if (!taxa_are_rows(physeq)) {
        otu_table(physeq) <-
          otu_table(
            t(as.matrix(unclass(
              physeq@otu_table
            ))),
            taxa_are_rows = TRUE
          )
      }
      df_physeq_interm <- cbind(
        physeq@otu_table,
        physeq@tax_table,
      )
      colnames(df_physeq_interm) <-
        c(
          sample_names(physeq),
          colnames(physeq@tax_table),
          "Reference Sequences"
        )

      df_physeq_interm <- as.data.frame(df_physeq_interm)

      if (write_sam_data) {
        sam_data <- data.frame(t(data.frame(unclass(
          physeq@sam_data
        ))))
        colnames(sam_data) <- sample_names(physeq)
        if (sam_data_first) {
          df_physeq <- dplyr::full_join(sam_data, df_physeq_interm)
          rownames(df_physeq) <-
            c(rownames(sam_data), rownames(df_physeq_interm))
        } else {
          df_physeq <- dplyr::full_join(df_physeq_interm, sam_data)
          rownames(df_physeq) <-
            c(rownames(df_physeq_interm), rownames(sam_data))
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
#' A wrapper of write_pq to save in all three possible formats
#'
#' @details
#' `r lifecycle::badge("maturing")`
#'
#' Write :
#' - 4 separate tables
#' - 1 table version
#' - 1 RData file
#'
#'
#' @inheritParams clean_pq
#' @param path a path to the folder to save the phyloseq object
#' @param ... Other arguments passed on to [write_pq()] or [utils::write.table()] function.
#' @return Build a folder (in path) with four csv tables (`refseq.csv`, `otu_table.csv`, `tax_table.csv`, `sam_data.csv`) + one
#'   table with all tables together + a rdata file (`physeq.RData`) that can be loaded using
#'   [base::load()] function + if present a phylogenetic tree in Newick format (`phy_tree.txt`)
#' @export
#' @author Adrien Taudière
#' @examples
#' \dontrun{
#' data(data_fungi)
#' save_pq(data_fungi, path = "phyloseq")
#' }
#' @seealso [MiscMetabar::write_pq()]
save_pq <- function(physeq, path = NULL, ...) {
  write_pq(physeq,
    path = path,
    rdata = TRUE,
    one_file = TRUE,
    ...
  )
  write_pq(physeq,
    path = path,
    rdata = FALSE,
    one_file = FALSE,
    ...
  )
}

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
#' @param sep_csv (default tabulation) separator for column
#' @param ... Other arguments passed on to [utils::write.table()] function.
#' @return One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#' and if present a phy_tree in Newick format. At least the otu_table.csv need to be present.
#' @export
#'
#' @examples
#' \dontrun{
#' read_pq(path = "phyloseq_data")
#' read_pq(path = "phyloseq_data", taxa_are_rows = TRUE)
#' }
#'
read_pq <- function(
        path = NULL,
        taxa_are_rows = FALSE,
        sam_names = NULL,
        sep_csv = "\t",
        ...) {
  if (file.exists(paste0(path, "/otu_table.csv"))) {
    if (taxa_are_rows) {
      otu_table_csv <-
        as.matrix(utils::read.table(paste0(path, "/otu_table.csv"), sep = sep_csv))
      samp_names <- colnames(otu_table_csv)
      otu_table_csv <- apply(otu_table_csv, 2, as.numeric)
      table_otu <- otu_table(otu_table_csv, taxa_are_rows = TRUE)
      sample_names(table_otu) <- samp_names
      physeq <- phyloseq(table_otu)
    } else {
      otu_table_csv <-
        as.matrix(utils::read.table(paste0(path, "/otu_table.csv"), sep = sep_csv))
      samp_names <- rownames(otu_table_csv)
      otu_table_csv <- apply(otu_table_csv, 2, as.numeric)
      rownames(otu_table_csv) <- samp_names
      physeq <-
        phyloseq(otu_table(otu_table_csv, taxa_are_rows = FALSE))
    }
  }
  if (file.exists(paste0(path, "/refseq.csv"))) {
    dna <-
      Biostrings::DNAStringSet(utils::read.table(
        paste0(path, "/refseq.csv"),
        sep = sep_csv,
        row.names = NULL
      )[, 2])
    names(dna) <-
      utils::read.table(paste0(path, "/refseq.csv"),
        sep = sep_csv,
        row.names = NULL
      )[, 1]
    physeq <- phyloseq::merge_phyloseq(physeq, refseq(dna))
  }
  if (file.exists(paste0(path, "/tax_table.csv"))) {
    tax_table_csv <-
      utils::read.table(paste0(path, "/tax_table.csv"), sep = sep_csv)
    rownames(tax_table_csv) <- tax_table_csv[, 1]
    tax_table_csv <- as.matrix(tax_table_csv[, -1])
    physeq <-
      phyloseq::merge_phyloseq(physeq, tax_table(tax_table_csv))
  }
  if (file.exists(paste0(path, "/sam_data.csv"))) {
    sam_data_csv <-
      utils::read.table(paste0(path, "/sam_data.csv"), sep = sep_csv)
    rownames(sam_data_csv) <- sam_data_csv[, 1]
    physeq <-
      phyloseq::merge_phyloseq(physeq, sample_data(sam_data_csv))
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
#' @description
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
#' @param  keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#' @return a list of for object
#' - "new_physeq": The new phyloseq object (class physeq)
#' - "discrepancy_vector": A vector of discrepancy showing for each taxonomic
#'   level the proportion of identic value before and after lulu reclustering.
#'   A value of 0.6 stands for 60% of ASV before re-clustering have
#'   identical value after re-clustering. In other word, 40% of ASV are assigned
#'   to a different taxonomic
#'   value. NA value are not counted as discrepancy.
#' - "res_lulu": A list of the result from the lulu function
#' - "merged_ASV": the data.frame used to merged ASV
#'
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
                    clean_pq = FALSE,
                    keep_temporary_files = FALSE) {
  verify_pq(physeq)
  if (is.null(physeq@refseq)) {
    stop("The phyloseq object do not contain a @refseq slot")
  }
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

  message("Lulu algorithm")
  res_lulu <-
    lulu(data.frame(t(physeq@otu_table)), match_list)

  if (!keep_temporary_files) {
    if (file.exists("temp.fasta")) {
      unlink("temp.fasta")
    }
    if (file.exists("cluster.fasta")) {
      unlink("cluster.fasta")
    }
    if (file.exists("temp.uc")) {
      unlink("temp.uc")
    }
    if (file.exists("match_list.txt")) {
      unlink("match_list.txt")
    }
  }
  merged <- res_lulu$otu_map[res_lulu$otu_map$curated == "merged", ]
  merged <- merged[rownames(merged) != merged$parent_id, ]

  test_vector <- vector(mode = "logical")
  for (tax_rank in colnames(physeq@tax_table)) {
    test <-
      physeq@tax_table[rownames(merged), tax_rank] == physeq@tax_table[merged$parent_id, tax_rank]
    test_vector <-
      c(
        test_vector,
        sum(test, na.rm = TRUE) / length(stats::na.exclude(test))
      )
  }

  names(test_vector) <- colnames(physeq@tax_table)

  new_physeq <-
    prune_taxa(
      taxa_names(physeq) %in% rownames(res_lulu$curated_table),
      physeq
    )
  new_physeq@otu_table <-
    otu_table(t(res_lulu$curated_table), taxa_are_rows = FALSE)
  sample_names(new_physeq) <- sample_names(physeq)

  if (verbose) {
    message(paste(
      "The number of taxa decrease from ",
      ntaxa(physeq),
      " to ",
      ntaxa(new_physeq),
      ".",
      sep = ""
    ))
    message(
      "See the discrepancy_vector to verify the degree of discrepancy in taxonomy due to lulu re-clustering."
    )
  }
  return(
    list(
      "new_physeq" = new_physeq,
      "discrepancy_vector" = test_vector,
      "res_lulu" = res_lulu,
      "merged_ASV" = merged
    )
  )
}
################################################################################


################################################################################
#' MUMU reclustering of class `physeq`
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' See https://www.nature.com/articles/s41467-017-01312-x for more information
#'  on the original method LULU. This is a wrapper of
#'  [mumu](https://github.com/frederic-mahe/mumu) a C++ re-implementation
#'  of LULU by Frédéric Mahé
#'
#' @inheritParams clean_pq
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param id (default: 0.84) id for --usearch_global.
#' @param vsearchpath path to vsearch.
#' @param mumupath path to mumu. See [mumu](https://github.com/frederic-mahe/mumu)
#'   for installation instruction
#' @param verbose (logical) if true, print some additional messages.
#' @param clean_pq (logical) if true, empty samples and empty ASV are discarded
#'   before clustering.
#' @param  keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#' @return a list of for object
#' - "new_physeq": The new phyloseq object (class physeq)
#' - "mumu_results": The log file of the mumu software. Run `man mumu` into
#'   bash to obtain details about columns' signification.
#'
#' @export
#' @examples
#' \dontrun{
#' data(data_fungi_sp_known)
#' mumu_pq(data_fungi_sp_known)
#' }
#' @author Frédéric Mahé
#'   & Adrien Taudière \email{adrien.taudiere@@zaclys.net}
#' @references
#' - MUMU: \url{https://github.com/frederic-mahe/mumu}
#' - VSEARCH can be downloaded from
#'  \url{https://github.com/torognes/vsearch}.
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [mumu](https://github.com/frederic-mahe/mumu/blob/main/CITATION.cff) and
#'   [lulu](https://www.nature.com/articles/s41467-017-01312-x) if you use this function
#'   for your work.
#'
mumu_pq <- function(
        physeq,
        nproc = 1,
        id = 0.84,
        vsearchpath = "vsearch",
        mumupath = "mumu",
        verbose = FALSE,
        clean_pq = TRUE,
        keep_temporary_files = FALSE) {
  verify_pq(physeq)
  if (is.null(physeq@refseq)) {
    stop("The phyloseq object do not contain a @refseq slot")
  }
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
      " -userfields query+target+id --maxaccepts 0 --query_cov 0.9 --maxhits 10",
      " -id ",
      id,
      "  --userout match_list.txt"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  otu_tab <-
    data.frame(unclass(clean_pq(physeq, force_taxa_as_rows = TRUE)@otu_table))
  otu_tab <- cbind("ASV" = rownames(otu_tab), otu_tab)
  write.table(
    otu_tab,
    "otu_table.csv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  message("Mumu algorithm")
  system2(
    mumupath,
    paste(
      "--otu_table otu_table.csv",
      "--match_list match_list.txt",
      "--log log.txt",
      "--new_otu_table new_OTU.tablemumu"
    )
  )

  res_mumu <- read.delim("new_OTU.tablemumu")
  new_otu_tab <- otu_table(res_mumu[, -1], taxa_are_rows = TRUE)
  taxa_names(new_otu_tab) <- res_mumu[, 1]

  new_physeq <-
    prune_taxa(taxa_names(physeq) %in% taxa_names(new_otu_tab), physeq)
  new_physeq@otu_table <-
    otu_table(t(new_otu_tab), taxa_are_rows = FALSE)
  if (nsamples(new_physeq) != nsamples(physeq)) {
    stop(
      "There is a different number of samples before and after mumu algorithm.
         This may be due to empty samples. You may try to rerun mumu_pq()
         using clean_pq = TRUE."
    )
  }
  sample_names(new_physeq) <- sample_names(physeq)

  if (verbose) {
    message(paste(
      "The number of taxa decrease from ",
      ntaxa(physeq),
      " to ",
      ntaxa(new_physeq),
      ".",
      sep = ""
    ))
    message(
      "See the log slot to verify the degree of discrepancy in taxonomy due to mumu re-clustering."
    )
  }

  result_mumu <- read.delim("log.txt")
  colnames(result_mumu) <-
    c(
      "Query_ASV",
      "Potential parent",
      "Similarity_percent",
      "Ab_query",
      "Ab_parent",
      "Overlap_ab_query",
      "Overlap_ab_parent",
      "Incidence_query",
      "Incidence_parent",
      "Smallest_ab_ratio",
      "Sum_ab_ratio",
      "Average_ratio",
      "Average_non_null_ratio",
      "Relative_cooccurence_value",
      "Status"
    )

  if (!keep_temporary_files) {
    if (file.exists("temp.fasta")) {
      unlink("temp.fasta")
    }
    if (file.exists("cluster.fasta")) {
      unlink("cluster.fasta")
    }
    if (file.exists("temp.uc")) {
      unlink("temp.uc")
    }

    if (file.exists("log.txt")) {
      unlink("temp.uc")
    }
    if (file.exists("match_list.txt")) {
      unlink("match_list.txt")
    }
    if (file.exists("otu_table.csv")) {
      unlink("otu_table.csv")
    }
    if (file.exists("new_OTU.tablemumu")) {
      unlink("new_OTU.tablemumu")
    }
  }

  return(list(
    "new_physeq" = new_physeq,
    "mumu_results" = result_mumu
  ))
}
################################################################################




################################################################################
#' Verify the validity of a phyloseq object
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' Mostly for internal use in MiscMetabar functions.
#'
#' @inheritParams clean_pq
#' @param verbose (logical, default FALSE) If true, prompt some warnings.
#' @param  min_nb_seq_sample (numeric) Only used if verbose = TRUE.
#'   Minimum number of sequences per samples to not show warning.
#' @param  min_nb_seq_taxa (numeric) Only used if verbose = TRUE.
#'   Minimum number of sequences per taxa to not show warning.
#' @return Nothing if the phyloseq object is valid. An error in the other case.
#'  Warnings if verbose = TRUE
#' @export
#'
verify_pq <- function(
    physeq,
    verbose = FALSE,
    min_nb_seq_sample = 500,
    min_nb_seq_taxa = 1) {
  if (!methods::validObject(physeq) ||
    !inherits(physeq, "phyloseq")) {
    stop("The physeq argument is not a valid phyloseq object.")
  }
  if (verbose) {
    if (min(sample_sums(physeq)) < min_nb_seq_sample) {
      warning(paste0(
        "At least one of your sample contains less than ",
        min_nb_seq_sample,
        " sequences."
      ))
    }
    if (min(sample_sums(physeq)) < min_nb_seq_sample) {
      warning(paste0(
        "At least one of your taxa is represent by less than ",
        min_nb_seq_taxa,
        " sequences."
      ))
    }
    if (sum(is.na(physeq@sam_data)) > 0) {
      warning("At least one of your samples metadata columns contains NA.")
    }
  }
}
################################################################################


################################################################################
#' Subset samples using a conditional boolean vector.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The main objective of this function is to complete the
#' [phyloseq::subset_samples()] function by propose a more easy
#' (but more prone to error) way of subset_samples.
#' It replace the subsetting expression which used the name of the variable
#' in the sam_data by a boolean vector.
#'
#' Warnings: you must verify the result of this function as the
#' boolean condition must match the order of samples in the `sam_data`
#' slot.
#'
#' This function is robust when you use the sam_data slot of the phyloseq object
#' used in physeq (see examples)
#'
#' @inheritParams clean_pq
#' @param condition A boolean vector to subset samples. Length must fit
#'   the number of samples
#'
#' @examples
#' data(data_fungi)
#' cond_samp <- grepl("A1", data_fungi@sam_data[["Sample_names"]])
#' subset_samples_pq(data_fungi, cond_samp)
#'
#' subset_samples_pq(data_fungi, data_fungi@sam_data[["Height"]] == "Low")
#'
#' @return a new phyloseq object
#' @export
#'
subset_samples_pq <- function(physeq, condition) {
  if (length(condition) != nsamples(physeq)) {
    stop("Length of condition is different from the number of samples.")
  }
  if (is.null(sample_data(physeq))) {
    message("Nothing subset. No sample_data in physeq.\n")
    return(physeq)
  } else {
    old_DF <- as(sample_data(physeq), "data.frame")
    new_DF <- old_DF[condition, ]
    if (inherits(physeq, "sample_data")) {
      return(sample_data(new_DF))
    } else {
      sample_data(physeq) <- sample_data(new_DF)
      return(physeq)
    }
  }
}
################################################################################

################################################################################
#' Subset taxa using a conditional named boolean vector.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The main objective of this function is to complete the
#' [phyloseq::subset_taxa()] function by propose a more easy way of
#' subset_taxa using a named boolean vector. Names must match taxa_names.
#'
#'
#' @inheritParams clean_pq
#' @param condition A named boolean vector to subset taxa. Length must fit
#'   the number of taxa and names must match taxa_names. Can also be a
#'   condition using a column of the tax_table slot (see examples). If
#'   the order of condition is the same as taxa_names(physeq),
#'   you can use the parameter `taxa_names_from_physeq = TRUE`.
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples are discarded after subsetting ASV
#' @param verbose (logical) Informations are printed
#' @param taxa_names_from_physeq (logical) If set to TRUE, rename the
#'   condition vector using taxa_names(physeq). Carefully check the result
#'   of this function if you use this parameter. No effect if the condition
#'   is of class `tax_table`.
#' @examples
#' data(data_fungi)
#' subset_taxa_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota")
#'
#' cond_taxa <- grepl("Endophyte", data_fungi@tax_table[, "Guild"])
#' names(cond_taxa) <- taxa_names(data_fungi)
#' subset_taxa_pq(data_fungi, cond_taxa)
#'
#' subset_taxa_pq(data_fungi, grepl("mycor", data_fungi@tax_table[, "Guild"]),
#'   taxa_names_from_physeq = TRUE
#' )
#'
#' @return a new phyloseq object
#' @export
#'
subset_taxa_pq <- function(
        physeq,
        condition,
        verbose = TRUE,
        clean_pq = TRUE,
        taxa_names_from_physeq = FALSE) {
  if (inherits(condition, "taxonomyTable")) {
    condition_temp <- as.vector(condition)
    names(condition_temp) <- rownames(condition)
    condition <- condition_temp
  } else {
    if (taxa_names_from_physeq) {
      names(condition) <- taxa_names(physeq)
    }
  }

  if (!sum(names(condition) %in% taxa_names(physeq)) == length(condition)) {
    stop(paste(
      "Some names in condition do not fit taxa_names of physeq : ",
      paste(names(condition)[!names(condition) %in% taxa_names(physeq)],
        collapse = "/"
      )
    ))
  }

  new_physeq <- physeq

  if (!taxa_are_rows(new_physeq)) {
    new_physeq@otu_table <-
      otu_table(t(new_physeq@otu_table), taxa_are_rows = TRUE)
    taxa_are_rows(new_physeq) <- TRUE
  }

  cond <- condition[match(taxa_names(new_physeq), names(condition))]
  cond[is.na(cond)] <- FALSE

  old_MA <-
    as(otu_table(new_physeq, taxa_are_rows = TRUE), "matrix")
  new_MA <- old_MA[cond, ]

  if (!is.matrix(new_MA)) {
    new_MA <- as.matrix(new_MA)
    new_otu_table <- otu_table(new_MA, taxa_are_rows = TRUE)
    sample_names(new_otu_table) <- sample_names(new_physeq)
  } else {
    new_otu_table <- otu_table(new_MA, taxa_are_rows = TRUE)
  }

  otu_table(new_physeq) <- new_otu_table

  if (clean_pq) {
    new_physeq <- clean_pq(new_physeq, verbose = TRUE)
  }

  if (verbose) {
    message(paste("Number of non-matching ASV", sum(is.na(
      match(taxa_names(physeq), names(condition))
    ))))
    message(paste("Number of matching ASV", sum(!is.na(
      match(taxa_names(physeq), names(condition))
    ))))
    message(paste(
      "Number of filtered-out ASV",
      ntaxa(physeq) - ntaxa(new_physeq)
    ))
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
#'
#' @author Adrien Taudière
#'
#' @examples
#' data(data_fungi)
#' A8_005 <- select_one_sample(data_fungi, "A8-005_S4_MERGED.fastq.gz")
#' A8_005
select_one_sample <- function(physeq, sam_name, silent = FALSE) {
  if (sum(sample_names(physeq) %in% sam_name) == 0) {
    stop(
      paste0(
        "The sample ",
        sam_name,
        " is not present in the names of samples of your phyloseq physeq object.
        You may use the sample_names() function."
      )
    )
  }
  cl_sam <-
    clean_pq(subset_samples_pq(physeq, sample_names(physeq) == sam_name),
      silent = TRUE
    )

  if (!silent) {
    message(
      paste0(
        "You select 1 of ",
        nsamples(physeq),
        " samples and conserved ",
        ntaxa(cl_sam),
        " out of ",
        ntaxa(physeq),
        " taxa represented by ",
        sum(cl_sam@otu_table),
        " sequences (out of ",
        sum(physeq@otu_table),
        " sequences [",
        perc(sum(cl_sam@otu_table), sum(physeq@otu_table)),
        "%])"
      )
    )
  }

  return(cl_sam)
}
################################################################################


################################################################################
#' Add new taxonomic rank to a phyloseq object.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' One of main use of this function is to add taxonomic assignment from
#' a new database.
#'
#' @inheritParams clean_pq
#' @param ref_fasta (required) A link to a database.
#'   Pass on to `dada2::assignTaxonomy`.
#' @param suffix (character) The suffix to name the new columns.
#'   If set to NULL (the default), the basename of the file reFasta
#'   is used.
#' @param ... Others arguments pass on to `dada2::assignTaxonomy`.
#' @return a physeq object with a larger slot tax_table
#'
#' @export
#' @examples
#' # example code
#'
#' @author Adrien Taudière
#'
add_new_taxonomy_pq <- function(physeq, ref_fasta, suffix = NULL, ...) {
  if (is.null(suffix)) {
    suffix <- basename(ref_fasta)
  }
  tax_tab <-
    dada2::assignTaxonomy(physeq@refseq, refFasta = ref_fasta, ...)
  colnames(tax_tab) <-
    make.unique(paste0(colnames(tax_tab), "_", suffix))
  new_tax_tab <- tax_table(cbind(physeq@tax_table, tax_tab))
  new_physeq <- physeq
  tax_table(new_physeq) <- new_tax_tab
  return(new_physeq)
}
################################################################################


################################################################################
#' Summarize information from sample data in a table
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [gtsummary::tbl_summary()] function in the case of `physeq`
#'   object.
#'
#' @inheritParams clean_pq
#' @param remove_col_unique_value (logical, default TRUE) Do we remove
#'  informative columns (categorical column with one value per samples),
#'   e.g. samples names ?
#' @param ... Others arguments pass on to [gtsummary::tbl_summary()].
#' @return a physeq object with a larger slot tax_table
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' tbl_sum_samdata(data_fungi) %>%
#'   gtsummary::as_kable()
#'
#' summary_samdata <- tbl_sum_samdata(data_fungi,
#'   include = c("Time", "Height"),
#'   type = list(Time ~ "continuous2", Height ~ "categorical"),
#'   statistic = list(Time ~ c("{median} ({p25}, {p75})", "{min}, {max}"))
#' )
#'
#' data(enterotype)
#'
#' summary_samdata <- tbl_sum_samdata(enterotype)
#' summary_samdata <- tbl_sum_samdata(enterotype, include = !contains("SampleId"))
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `gtsummary::tbl_summary()` if you
#'   use this function.

tbl_sum_samdata <- function(physeq, remove_col_unique_value = TRUE, ...) {
  tbl <- tibble(data.frame(physeq@sam_data))
  if (remove_col_unique_value) {
    tbl <- tbl[, !apply(tbl, 2, function(x) {
      length(unique(x)) == nrow(tbl) && is.character(x)
    })]
  }
  tbl_sum <- tbl %>% gtsummary::tbl_summary(...)
  return(tbl_sum)
}
################################################################################
#' Add information about Guild for FUNGI the FUNGuild databse
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Please cite this [publication](https://doi.org/10.1016/j.funeco.2015.06.006).
#'
#' @inheritParams clean_pq
#' @param taxLevels Name of the 7 columns in tax_table required by funguild
#'
#' @return A new object of class `physeq` with Guild information added to
#'   `tax_table` slot
#' @export
#' @author Adrien Taudière
#' @examples
#' data(data_fungi)
#' df <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 5000)
#' \dontrun{
#' df <- add_funguild_info(df,
#'   taxLevels = c(
#'     "Domain", "Phylum", "Class", "Order",
#'     "Family", "Genus", "Species"
#'   )
#' )
#' sort(table(df@tax_table[, "guild"]), decreasing = TRUE)
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `FUNGuildR` package and the associate
#'   [publication](https://www.sciencedirect.com/science/article/abs/pii/S1754504815000847) if you
#'   use this function.
#' @seealso [plot_guild_pq()]

add_funguild_info <- function(
        physeq,
        taxLevels = c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species"
        )) {
  tax_tab <- physeq@tax_table
  FUNGuild_assign <-
    funguild_assign(data.frame(
      "Taxonomy" =
        apply(tax_tab[, taxLevels], 1,
          paste,
          collapse = ";"
        )
    ))
  tax_tab <-
    as.matrix(cbind(
      tax_tab,
      FUNGuild_assign
    ))
  physeq@tax_table <- tax_table(tax_tab)
  return(physeq)
}




################################################################################
#' Plot information about Guild from tax_table slot previously
#' created with [add_funguild_info()]
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams clean_pq
#' @param levels_order (Default NULL) A character vector to
#'   reorder the levels of guild. See examples.
#' @param clean_pq (logical, default TRUE): Does the phyloseq
#'   object is cleaned using the [clean_pq()] function?
#' @param ... other params for be passed on to
#'   [clean_pq()] function
#' @return A ggplot2 object
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' \dontrun{
#' data(data_fungi)
#' df <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 5000)
#' df <- add_funguild_info(df,
#'   taxLevels = c(
#'     "Domain", "Phylum", "Class", "Order",
#'     "Family", "Genus", "Species"
#'   )
#' )
#' p <- plot_guild_pq(df)
#' library("patchwork")
#' (plot_guild_pq(subset_samples(df, Height == "Low"),
#'   levels_order = p$data$Guild[order(p$data$nb_seq)]
#' ) + theme(legend.position = "none")) +
#'   (plot_guild_pq(subset_samples(df, Height == "High"),
#'     levels_order = p$data$Guild[order(p$data$nb_seq)]
#'   ) + ylab("") + theme(axis.text.y = element_blank()))
#' }
#' @seealso [add_funguild_info()]

plot_guild_pq <-
  function(
        physeq,
        levels_order = NULL,
        clean_pq = TRUE,
        ...) {
    if (clean_pq) {
      physeq <- clean_pq(physeq, ...)
    }
    guilds <-
      data.frame(sort(table(strsplit(
        paste(
          physeq@tax_table[, "guild"]
          [physeq@tax_table[, "confidenceRanking"] %in%
              c("Highly Probable", "Probable")],
          collapse = "-"
        ),
        split = "-"
      ))))

    guilds$Var1 <- as.vector(guilds$Var1)
    guilds <- guilds[guilds$Var1 != "NA", ]
    guilds <- guilds[guilds$Var1 != "NULL", ]
    guilds <- guilds[guilds$Var1 != "", ]

    # Number of sequences per guild
    nb_seq_by_guild <- c()
    for (i in seq(1, length(guilds$Var1))) {
      nb_seq_by_guild[i] <-
        sum(taxa_sums(physeq@otu_table)[grepl(
          guilds$Var1[i],
          physeq@tax_table[, "guild"]
        )])
    }
    names(nb_seq_by_guild) <- guilds$Var1
    guilds$seq <- nb_seq_by_guild

    names(guilds) <- c("Guild", "nb_asv", "nb_seq")
    guilds$nb_seq <- as.numeric(guilds$nb_seq)
    guilds$nb_asv <- as.numeric(guilds$nb_asv)

    guilds$Guild <- factor(as.vector(guilds$Guild),
      levels = guilds$Guild[order(guilds$nb_seq)]
    )


    COLORS <- rep("Others", nrow(guilds))
    COLORS[grepl("Sapro", guilds$Guild)] <- "Sapro"
    COLORS[grepl("Parasite", guilds$Guild)] <- "Parasite/pathogen"
    COLORS[grepl("Pathog", guilds$Guild)] <- "Parasite/pathogen"
    COLORS[grepl("ycorrh", guilds$Guild)] <- "Mutualist"
    COLORS[grepl("Lichen", guilds$Guild)] <- "Mutualist"
    COLORS[grepl("Endophy", guilds$Guild)] <- "Mutualist"

    guilds$colors <- COLORS
    guilds <- rbind(
      guilds,
      data.frame(
        "Guild" = "All ASV",
        "nb_asv" = ntaxa(physeq),
        "nb_seq" = sum(physeq@otu_table),
        "colors" = "ALL"
      )
    )
    guilds <- guilds[order(guilds$nb_seq), ]
    if (!is.null(levels_order)) {
      guilds$Guild <- factor(guilds$Guild, levels = levels_order)
    }

    ggplot(
      guilds,
      aes(
        y = Guild,
        x = log10(nb_seq),
        fill = colors
      )
    ) +
      geom_bar(stat = "identity") +
      annotation_logticks(sides = "b", alpha = 0.5) +
      ylab("GUILD by FUNGuild") +
      scale_fill_manual("Guild",
        values = c(
          "gray", "Olivedrab", "cyan4", "tomato3",
          "lightpink4"
        )
      ) +
      geom_text(aes(label = nb_asv, x = log10(nb_seq) + 0.2),
        family = "serif"
      ) +
      geom_text(aes(label = nb_seq, x = log10(nb_seq) / 2),
        family = "mono",
        col = "white"
      )
  }

################################################################################

################################################################################
#' Build phylogenetic trees from refseq slot of a phyloseq object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function build tree phylogenetic tree and if nb_bootstrap is
#' set, it build also the 3 corresponding bootstrapped tree.
#'
#' Default parameters are based on https://doi.org/10.12688/f1000research.8986.2
#' and phangorn vignette [Estimating phylogenetic trees with phangorn](https://klausvigo.github.io/phangorn/articles/Trees.html). You should understand your data, especially the markers,
#' before using this function.
#'
#' Note that phylogenetic reconstruction with markers used for metabarcoding are
#' not robust. You must verify the robustness of your phylogenetic tree using
#' taxonomic classification (see vignette [Tree visualization](https://adrientaudiere.github.io/MiscMetabar/articles/tree_visualization.html)) and bootstrap or multi-tree visualization
#'
#' @inheritParams clean_pq
#' @param nb_bootstrap (default 0): If a positive number is set,
#'   the function also build 3 bootstrapped trees using `nb_bootstrap`
#'   bootstrap samples
#' @param model allows to choose an amino acid models or nucleotide model,
#'   see [phangorn::optim.pml()] for more details
#' @param optInv 	Logical value indicating whether topology gets optimized
#'  (NNI). See [phangorn::optim.pml()] for more details
#' @param optGamma Logical value indicating whether gamma rate parameter gets
#'  optimized. See [phangorn::optim.pml()] for more details
#' @param rearrangement type of tree tree rearrangements to perform, one of
#'  "NNI", "stochastic" or "ratchet"
#'   see [phangorn::optim.pml()] for more details
#' @param control A list of parameters for controlling the fitting process.
#'   see [phangorn::optim.pml()] for more details
#' @param optNni Logical value indicating whether topology gets optimized (NNI).
#'   see [phangorn::optim.pml()] for more details
#' @param multicore	(logical) whether models should estimated in parallel.
#'   see [phangorn::bootstrap.pml()] for more details
#' @param ... other params for be passed on to
#'   [phangorn::optim.pml()] function
#'
#' @return A list of phylogenetic tree
#' @export
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `phangorn` package if you
#'   use this function.
#' @examples
#' library("phangorn")
#' data(data_fungi)
#' df <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 9000)
#' df_tree <- build_phytree_pq(df, nb_bootstrap = 5)
#' plot(df_tree$UPGMA)
#' plotBS(df_tree$UPGMA, df_tree$UPGMA_bs, main = "UPGMA")
#' plot(df_tree$NJ, "unrooted")
#' plot(df_tree$ML)
#' # plotBS(df_tree$ML_bs)
#' plotBS(df_tree$ML$tree, df_tree$ML_bs, p = 20, frame = "circle")
#' plotBS(
#'   df_tree$ML$tree,
#'   df_tree$ML_bs,
#'   p = 20,
#'   frame = "circle",
#'   method = "TBE"
#' )
#' plot(consensusNet(df_tree$ML_bs))
#' plot(consensusNet(df_tree$NJ_bs))
#' ps_tree <- merge_phyloseq(df, df_tree$ML$tree)
build_phytree_pq <- function(
        physeq,
        nb_bootstrap = 0,
        model = "GTR",
        optInv = TRUE,
        optGamma = TRUE,
        rearrangement = "NNI",
        control = phangorn::pml.control(trace = 0),
        optNni = TRUE,
        multicore = FALSE,
        ...) {
  seqs <- physeq@refseq
  alignment <-
    DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor = NA)

  phang.align <-
    phangorn::phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- phangorn::dist.ml(phang.align)
  treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
  treeUPGMA <-
    phangorn::upgma(dm) # Note, tip order != sequence order
  fit <- phangorn::pml(treeNJ, data = phang.align)
  ## negative edges length changed to 0!
  fitGTR <- update(fit, k = 4, inv = 0.2)
  tree_ML <-
    phangorn::optim.pml(
      fitGTR,
      model = model,
      optInv = optInv,
      optGamma = optGamma,
      rearrangement = rearrangement,
      control = control,
      optNni = optNni,
      ...
    )
  if (nb_bootstrap > 0) {
    treeUPGMA_bs <-
      phangorn::bootstrap.phyDat(phang.align,
        function(x) {
          phangorn::upgma(phangorn::dist.ml(x))
        },
        bs = nb_bootstrap
      )
    if (rearrangement == "NNI") {
      tree_ML_bs <- phangorn::bootstrap.pml(
        tree_ML,
        bs = nb_bootstrap,
        multicore = multicore,
        rearrangement = "NNI",
        ...
      )
    } else if (rearrangement == "stochastic") {
      tree_ML_bs <- phangorn::bootstrap.pml(
        tree_ML,
        bs = nb_bootstrap,
        multicore = multicore,
        rearrangement = "stochastic",
        ...
      )
    } else if (rearrangement == "ratchet") {
      tree_ML_bs <- phangorn::bootstrap.pml(
        tree_ML,
        bs = nb_bootstrap,
        multicore = multicore,
        rearrangement = "ratchet",
        ...
      )
    } else {
      stop("rearrangement parameter one of the three value 'stochastic',
       'NNI' or 'ratchet'")
    }
    treeNJ_bs <- phangorn::bootstrap.phyDat(phang.align,
      function(x) {
        phangorn::NJ(phangorn::dist.ml(x))
      },
      bs = nb_bootstrap
    )
    return(
      list(
        "UPGMA" = treeUPGMA,
        "NJ" = treeNJ,
        "ML" = tree_ML,
        "UPGMA_bs" = treeUPGMA_bs,
        "NJ_bs" = treeNJ_bs,
        "ML_bs" = tree_ML_bs
      )
    )
  } else {
    return(list(
      "UPGMA" = treeUPGMA,
      "NJ" = treeNJ,
      "ML" = tree_ML
    ))
  }
}
################################################################################

################################################################################
#' Test if the mean number of sequences by samples is link to the modality of
#' a factor
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' The aim of this function is to provide a warnings if samples depth significantly
#' vary among the modalities of a factor present in the `sam_data` slot.
#'
#' This function apply a Kruskal-Wallis rank sum test to the number of sequences
#' per samples in function of the factor `fact`.
#'
#' @inheritParams clean_pq
#' @param fact (required): Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}.
#' @param boxplot (logical) Do you want to plot boxplot?
#'
#' @return The result of a Kruskal-Wallis rank sum test
#' @export
#' @author Adrien Taudière
#' @importFrom stats kruskal.test
#' @examples
#' data(data_fungi)
#' are_modality_even_depth(data_fungi, "Time")$p.value
#' are_modality_even_depth(rarefy_even_depth(data_fungi), "Time")$p.value
#' are_modality_even_depth(data_fungi, "Height", boxplot = TRUE)
are_modality_even_depth <- function(physeq, fact, boxplot = FALSE) {
  nb_seq <- sample_sums(physeq)
  fact <- factor(unclass(physeq@sam_data[, fact])[[1]])
  res <- kruskal.test(nb_seq ~ fact)
  if (boxplot) {
    boxplot(nb_seq ~ fact)
  }
  if (length(unique(tapply(nb_seq, fact, mean))) == 1) {
    message("All modality were undoubtedly rarefy in the physeq object.")
    res$p.value <- 1
  }
  return(res)
}
################################################################################


################################################################################
#' Reorder taxa in otu_table/tax_table/refseq slot of a phyloseq object
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Note that the taxa order in a physeq object with a tree is locked by
#' the order of leaf in the phylogenetic tree.
#'
#' @inheritParams clean_pq
#' @param names_ordered (required): Names of the taxa (must be the same
#'   as taxa in `taxa_names(physeq)`) in a given order
#' @param remove_phy_tree (logical, default FALSE) If TRUE, the phylogenetic
#'   tree is removed. It is
#' @return A phyloseq object
#' @export
#' @author Adrien Taudière
#' @examples
#' data(data_fungi)
#' data_fungi_ordered_by_genus <- reorder_taxa_pq(
#'   data_fungi,
#'   taxa_names(data_fungi)[order(as.vector(data_fungi@tax_table[, "Genus"]))]
#' )
#'
#' data_fungi_asc_ordered_by_abundance <- reorder_taxa_pq(
#'   data_fungi,
#'   taxa_names(data_fungi)[order(taxa_sums(data_fungi))]
#' )
reorder_taxa_pq <- function(physeq, names_ordered, remove_phy_tree = FALSE) {
  new_physeq <- physeq

  if (!is.null(phy_tree(new_physeq, FALSE))) {
    if (remove_phy_tree) {
      message("Removing phylogenetic tree!")
      new_physeq@phy_tree <- NULL
    } else {
      stop("The taxa order in a physeq object with a tree is locked by
      the order of leaf in the phylogenetic tree. You could use args
      remove_phy_tree = TRUE.")
    }
  }

  if (sum(!names_ordered %in% taxa_names(new_physeq)) > 0) {
    stop(
      "The taxa names in physeq and in the args names_ordered are not the same.
      You can try identify them usign match(names_ordered, taxa_names(physeq))"
    )
  }
  order_taxa_names <- match(names_ordered, taxa_names(new_physeq))
  if (!is.null(otu_table(new_physeq, FALSE))) {
    if (taxa_are_rows(new_physeq)) {
      new_physeq@otu_table <-
        otu_table(new_physeq, taxa_are_rows = TRUE)[order_taxa_names, ]
    } else {
      new_physeq@otu_table <-
        otu_table(new_physeq, taxa_are_rows = FALSE)[, order_taxa_names]
    }
  }
  if (!is.null(tax_table(new_physeq, FALSE))) {
    new_physeq@tax_table <-
      tax_table(new_physeq)[order_taxa_names, ]
  }
  if (!is.null(refseq(new_physeq, FALSE))) {
    new_physeq@refseq <-
      refseq(new_physeq)[order_taxa_names]
  }
  return(new_physeq)
}
################################################################################
