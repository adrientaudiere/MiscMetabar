if (getRversion() >= "2.15.1") {
  utils::globalVariables(".")
}

################################################################################
#' Add dna in `refseq` slot of a `physeq` object using taxa names and renames taxa
#'   using prefix_taxa_names and number (default Taxa_1, Taxa_2 ...)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' Useful in targets bioinformatic pipeline.
#'
#' @inheritParams clean_pq
#' @param prefix_taxa_names (default "Taxa_"): the prefix of taxa names (eg. "ASV_" or "OTU_")
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with `refseq` slot and new
#'   taxa names
#' @author Adrien Taudière
#' @export

add_dna_to_phyloseq <- function(physeq, prefix_taxa_names = "Taxa_") {
  verify_pq(physeq)
  dna <- Biostrings::DNAStringSet(phyloseq::taxa_names(physeq))
  names(dna) <- phyloseq::taxa_names(physeq)
  physeq <- phyloseq::merge_phyloseq(physeq, dna)
  phyloseq::taxa_names(physeq) <-
    paste0(prefix_taxa_names, seq(phyloseq::ntaxa(physeq)))
  return(physeq)
}
################################################################################



################################################################################
#'  Clean phyloseq object by removing empty samples and taxa
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'  In addition, this function check for discrepancy (and rename) between
#' (i) taxa names in refseq, taxonomy table and otu_table and between
#' (ii) sample names in sam_data and otu_table.
#'
#' @param physeq (required) a \code{\link[phyloseq]{phyloseq-class}} object obtained
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
#' @param reorder_taxa (logical) if TRUE the otu_table is ordered by the number of
#'   sequences of taxa (ASV, OTU) in descending order. Default to FALSE.
#' @param rename_taxa (logical) if TRUE, taxa (ASV, OTU) are renamed by their position
#'   in the OTU_table and prefix_taxa_names param (by default: Taxa_1, Taxa_2, ...).
#'   Default to FALSE. If rename taxa (ASV, OTU) is true,
#'   the taxa (ASV, OTU) names in verbose information can be misleading.
#' @param simplify_taxo (logical) if TRUE, correct the taxonomy_table using the
#'   `MiscMetabar::simplify_taxo()` function
#' @param prefix_taxa_names (default "Taxa_"): the prefix of taxa names (eg. "ASV_" or "OTU_")
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object
#' @export
#' @author Adrien Taudière
clean_pq <- function(physeq,
                     remove_empty_samples = TRUE,
                     remove_empty_taxa = TRUE,
                     clean_samples_names = TRUE,
                     silent = FALSE,
                     verbose = FALSE,
                     force_taxa_as_columns = FALSE,
                     force_taxa_as_rows = FALSE,
                     reorder_taxa = FALSE,
                     rename_taxa = FALSE,
                     simplify_taxo = FALSE,
                     prefix_taxa_names = "_Taxa") {
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

  if (reorder_taxa) {
    physeq <- reorder_taxa_pq(physeq, taxa_names(physeq)[order(taxa_sums(physeq), decreasing = TRUE)])
  }

  if (rename_taxa) {
    taxa_names(physeq) <- paste0(prefix_taxa_names, seq(1, ntaxa(physeq)))
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
      otu_table(t(as.matrix(unclass(
        physeq@otu_table
      ))), taxa_are_rows = FALSE)
    message("Taxa are now in columns.")
  }

  if (force_taxa_as_rows && !taxa_are_rows(physeq)) {
    otu_table(physeq) <-
      otu_table(t(as.matrix(unclass(
        physeq@otu_table
      ))), taxa_are_rows = TRUE)
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
    if (ntaxa(physeq) - ntaxa(new_physeq) != 0 ||
      nsamples(physeq) - nsamples(new_physeq) != 0) {
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
  }

  verify_pq(new_physeq)
  return(new_physeq)
}






################################################################################
#' Track the number of reads (= sequences), samples and cluster (e.g. ASV)
#' from various objects including dada-class and derep-class.
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
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
#' @param verbose (logical) If true, print some additional messages.
#' @param ... Additional arguments passed on to [clean_pq()] function.
#'
#' @return The number of sequences, clusters (e.g. OTUs, ASVs) and samples for
#'   each object.
#' @export
#' @author Adrien Taudière
#' @seealso [track_wkflow_samples()]
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' data(enterotype)
#' if (requireNamespace("pbapply")) {
#'   track_wkflow(list(data_fungi, enterotype), taxonomy_rank = c(3, 5))
#'   track_wkflow(list(
#'     "data FUNGI" = data_fungi,
#'     "fastq files forward" =
#'       unlist(list_fastq_files(system.file("extdata", package = "MiscMetabar"),
#'         paired_end = FALSE
#'       ))
#'   ))
#' }
track_wkflow <- function(list_of_objects,
                         obj_names = NULL,
                         clean_pq = FALSE,
                         taxonomy_rank = NULL,
                         verbose = TRUE,
                         ...) {
  if (verbose) {
    message("Compute the number of sequences")
  }
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
      if (verbose) {
        message(paste("Start object of class:", class(object), sep = " "))
      }
      if (inherits(object, "phyloseq")) {
        sum(object@otu_table)
      } else if (inherits(object, "matrix")) {
        sum(object, na.rm = TRUE)
      } else if (is.character(object[1]) &&
        length(object[1]) == 1 &&
        file.exists(object[1])) {
        if (summary(file(object[[1]]))$class == "gzfile") {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("zcat ", x, " | grep -c '^+$'", sep = ""), intern = TRUE))
          })
        } else if (grepl("\\.fastq$", object[1])) {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("cat ", x, " | grep -c '^+$'", sep = ""), intern = TRUE))
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
      if (verbose) {
        message(paste("Start object of class:", class(object), sep = " "))
      }
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
      if (verbose) {
        message(paste("Start object of class:", class(object), sep = " "))
      }
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
        if (verbose) {
          message(paste("Start object of class:", class(object), sep = " "))
        }
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
        if (verbose) {
          message(paste("Start object of class:", class(object), sep = " "))
        }
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

    rownames(track) <- c("nb_sequences", "nb_clusters", "nb_samples")
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
#' for each sample
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Contrary to [track_wkflow()], only phyloseq object are possible.
#' More information are available in the manual of the function [track_wkflow()]
#'
#' @param list_pq_obj (required) a list of object passed on to [track_wkflow()]
#'   Only phyloseq object will return value because information of sample is needed
#' @param ... Other args passed on to [track_wkflow()]
#'
#' @return A list of dataframe. cf [track_wkflow()] for more information
#'
#' @export
#' @author Adrien Taudière
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' tree_A10_005 <- subset_samples(data_fungi, Tree_name == "A10-005")
#' if (requireNamespace("pbapply")) {
#'   track_wkflow_samples(tree_A10_005)
#' }
track_wkflow_samples <- function(list_pq_obj, ...) {
  if (!inherits(list_pq_obj, "list")) {
    list_pq_obj <- list(list_pq_obj)
  }
  if (sum(!unlist(lapply(list_pq_obj, inherits, "phyloseq"))) != 0) {
    stop("At least one object in your list_pq_obj is not a phyloseq obj.")
  }
  sam_names <- unique(unlist(lapply(list_pq_obj, sample_names)))
  res <- vector("list", length(sam_names))
  names(res) <- sam_names
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
#'   or a list of DNA sequences
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' This function use the `merge_taxa_vec` function to  merge taxa into clusters.
#'
#' @inheritParams clean_pq
#' @param dna_seq You may directly use a character vector of DNA sequences
#'   in place of physeq args. When physeq is set, dna sequences take the value of
#'   `physeq@refseq`
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param method (default: clusterize)
#'   Set the clustering method.
#'   - `clusterize` use the [DECIPHER::Clusterize()] fonction,
#'   - `vsearch` use the vsearch software (https://github.com/torognes/vsearch)
#'     with arguments `--cluster_size` by default (see args `vsearch_cluster_method`)
#'     and `-strand both` (see args `vsearch_args`)
#'   - `swarm` use the swarm
#' @param vsearchpath (default: vsearch) path to vsearch
#' @param id (default: 0.97) level of identity to cluster
#' @param tax_adjust (Default 0) See the man page
#'   of [merge_taxa_vec()] for more details.
#'   To conserved the taxonomic rank of the most abundant taxa (ASV, OTU,...),
#'   set tax_adjust to 0 (default). For the moment only tax_adjust = 0 is
#'   robust
#' @param rank_propagation (logical, default FALSE). Do we propagate the
#' NA value from lower taxonomic rank to upper rank?
#' See the man page of [merge_taxa_vec()] for more details. 
#' @param vsearch_cluster_method (default: "--cluster_size) See other possible
#'   methods in the [vsearch manual](https://github.com/torognes/vsearch/) (e.g. `--cluster_size` or `--cluster_smallmem`)
#'   - `--cluster_fast` : Clusterize the fasta sequences in filename, automatically sort by decreasing sequence length beforehand.
#'   - `--cluster_size` : Clusterize the fasta sequences in filename, automatically sort by decreasing sequence abundance beforehand.
#'   - `--cluster_smallmem` : Clusterize the fasta sequences in filename without automatically modifying their order beforehand. Sequence are expected to be sorted by decreasing sequence length, unless *--usersort* is used.
#'     In that case you may set `vsearch_args` to vsearch_args = "--strand both --usersort"
#' @param vsearch_args (default : "--strand both") a one length character element defining other parameters to
#'   passed on to vsearch.
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - temp.fasta (refseq in fasta or dna_seq sequences)
#'   - cluster.fasta (centroid if method = "vsearch")
#'   - temp.uc (clusters if method = "vsearch")
#' @param swarmpath (default: swarm) path to swarm
#' @param d (default: 1) maximum number of differences allowed between two
#'   amplicons, meaning that two amplicons will be grouped if they have `d`
#'   (or less) differences
#' @param swarm_args (default : "--fastidious") a one length character
#'   element defining other parameters to  passed on to swarm See other possible
#'   methods in the [SWARM pdf manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)
#' @param method_clusterize (default "overlap") the method for the [DECIPHER::Clusterize()] method
#' @param ... Additional arguments passed on to [DECIPHER::Clusterize()]
#' @details This function use the `merge_taxa_vec` function to
#'   merge taxa into clusters. By default tax_adjust = 0. See the man page
#'   of [merge_taxa_vec()].
#'
#' @return A new object of class `physeq` or a list of cluster if dna_seq
#'   args was used.
#'
#' @examples
#' if (requireNamespace("DECIPHER")) {
#'   postcluster_pq(data_fungi_mini)
#' }
#' \donttest{
#' if (requireNamespace("DECIPHER")) {
#'   postcluster_pq(data_fungi_mini, method_clusterize = "longest")
#'
#'   if (MiscMetabar::is_swarm_installed()) {
#'     d_swarm <- postcluster_pq(data_fungi_mini, method = "swarm")
#'   }
#'   if (MiscMetabar::is_vsearch_installed()) {
#'     d_vs <- postcluster_pq(data_fungi_mini, method = "vsearch")
#'   }
#' }
#' }
#' @references
#'   VSEARCH can be downloaded from
#'   \url{https://github.com/torognes/vsearch}.
#'   More information in the associated publication
#'   \url{https://pubmed.ncbi.nlm.nih.gov/27781170}.
#' @seealso [vsearch_clustering()] and [swarm_clustering()]
#' @export
#' @author Adrien Taudière

postcluster_pq <- function(physeq = NULL,
                           dna_seq = NULL,
                           nproc = 1,
                           method = "clusterize",
                           id = 0.97,
                           vsearchpath = "vsearch",
                           tax_adjust = 0,
                           rank_propagation = FALSE,
                           vsearch_cluster_method = "--cluster_size",
                           vsearch_args = "--strand both",
                           keep_temporary_files = FALSE,
                           swarmpath = "swarm",
                           d = 1,
                           swarm_args = "--fastidious",
                           method_clusterize = "overlap",
                           ...) {
  if (inherits(physeq, "phyloseq")) {
    verify_pq(physeq)
    if (is.null(physeq@refseq)) {
      stop("The phyloseq object do not contain a @refseq slot")
    }
    dna <- Biostrings::DNAStringSet(physeq@refseq)
    if (!is.null(dna_seq)) {
      stop("You must use either physeq or dna_seq args but not both")
    }
  } else if (inherits(dna_seq, "character")) {
    dna <- Biostrings::DNAStringSet(dna_seq)
  } else {
    stop(
      "You must set the args physeq (object of class phyloseq) or
    dna_seq (character vector)."
    )
  }

  if (method == "clusterize") {
    ## Find clusters of ASVs to form the new OTUs
    clusters <- DECIPHER::Clusterize(
      dna,
      cutoff = 1 - id,
      # e.g. `cutoff = 0.03` for a 97% OTU
      processors = nproc,
      method = method_clusterize,
      ...
    )

    if (inherits(physeq, "phyloseq")) {
      new_obj <-
        merge_taxa_vec(physeq, clusters$cluster, tax_adjust = tax_adjust, rank_propagation = rank_propagation)
    } else if (inherits(dna_seq, "character")) {
      new_obj <- clusters
    } else {
      stop(
        "You must set the args physeq (object of class phyloseq) or
    dna_seq (character vector)."
      )
    }
  } else if (method == "vsearch") {
    new_obj <- vsearch_clustering(
      physeq = physeq,
      dna_seq = dna_seq,
      nproc = nproc,
      id = id,
      vsearchpath = vsearchpath,
      tax_adjust = tax_adjust,
      rank_propagation = rank_propagation,
      vsearch_cluster_method = vsearch_cluster_method,
      vsearch_args = vsearch_args,
      keep_temporary_files = keep_temporary_files
    )
  } else if (method == "swarm") {
    new_obj <- swarm_clustering(
      physeq = physeq,
      dna_seq = dna_seq,
      d = d,
      swarmpath = swarmpath,
      vsearch_path = vsearchpath,
      nproc = nproc,
      swarm_args = swarm_args,
      tax_adjust = tax_adjust,
      rank_propagation = rank_propagation,
      keep_temporary_files = keep_temporary_files
    )
  } else {
    stop("Method allows 2 values only : `clusterize`, `vsearch` or `swarm`")
  }
  return(new_obj)
}
################################################################################



################################################################################
#' @rdname postcluster_pq
#' @export
asv2otu <- postcluster_pq
################################################################################



################################################################################
#' Save phyloseq object in the form of multiple csv tables.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' This is the reverse function of [read_pq()].
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
#'   If set to TRUE, empty samples are discarded after subsetting taxa (ASV, OTU, ...)
#' @param reorder_taxa (logical) if TRUE the otu_table is ordered by the number of
#'   sequences of taxa (ASV, OTU, ...) (descending order). Default to TRUE. Only possible if clean_pq
#'   is set to TRUE.
#' @param rename_taxa reorder_taxa (logical) if TRUE, taxa (ASV, OTU, ...) are renamed by their position
#'   in the OTU_table (asv_1, asv_2, ...). Default to FALSE. Only possible if clean_pq
#'   is set to TRUE.
#' @param quote a logical value (default FALSE) or a numeric vector.
#'   If TRUE, any character or factor columns will be surrounded by
#'   double quotes.  If a numeric vector, its elements are taken
#'   as the indices of columns to quote.  In both cases, row and
#'   column names are quoted if they are written. If FALSE nothing is quoted.
#' @param sep_csv (default tabulation) separator for column
#' @param ... Additional arguments passed on to [utils::write.table()] function.
#' @return Build a folder (path) containing one to four csv tables
#'   (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#'   and if present a phy_tree in Newick format
#' @export
#' @author Adrien Taudière
#' @examples
#' write_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"))
#' write_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"), one_file = TRUE)
#' unlink(paste0(tempdir(), "/phyloseq"), recursive = TRUE)
#' @seealso [MiscMetabar::save_pq()]

write_pq <- function(physeq,
                     path = NULL,
                     rdata = FALSE,
                     one_file = FALSE,
                     write_sam_data = TRUE,
                     sam_data_first = FALSE,
                     clean_pq = TRUE,
                     reorder_taxa = FALSE,
                     rename_taxa = FALSE,
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
    reorder_taxa = reorder_taxa,
    rename_taxa = rename_taxa,
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
          otu_table(t(as.matrix(unclass(
            physeq@otu_table
          ))), taxa_are_rows = TRUE)
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
          otu_table(t(as.matrix(unclass(
            physeq@otu_table
          ))), taxa_are_rows = TRUE)
      }
      df_physeq_interm <- cbind(physeq@otu_table, physeq@tax_table, )
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Write :
#' - 4 separate tables
#' - 1 table version
#' - 1 RData file
#'
#' @inheritParams clean_pq
#' @param path a path to the folder to save the phyloseq object
#' @param ... Additional arguments passed on to [write_pq()] or [utils::write.table()] function.
#' @return Build a folder (in path) with four csv tables (`refseq.csv`, `otu_table.csv`, `tax_table.csv`, `sam_data.csv`) + one
#'   table with all tables together + a rdata file (`physeq.RData`) that can be loaded using
#'   [base::load()] function + if present a phylogenetic tree in Newick format (`phy_tree.txt`)
#' @export
#' @author Adrien Taudière
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' save_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"))
#' unlink(paste0(tempdir(), "/phyloseq"), recursive = TRUE)
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
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' This is the reverse function of [write_pq()].
#'
#' @param path (required) a path to the folder to read the phyloseq object
#' @param taxa_are_rows (default to FALSE) see ?phyloseq for details
#' @param sam_names The name of the variable (column) in sam_data.csv to rename
#'   samples. Note that if you use [write_phyloseq()] function to save your
#'   physeq object, you may use sam_names = "X" to rename the samples names
#'   as before.
#' @param sep_csv (default tabulation) separator for column
#' @param ... Additional arguments passed on to [utils::write.table()] function.
#' @return One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#' and if present a phy_tree in Newick format. At least the otu_table.csv need to be present.
#' @export
#' @author Adrien Taudière
#' @examples
#' write_pq(data_fungi, path = paste0(tempdir(), "/phyloseq"))
#' read_pq(path = paste0(tempdir(), "/phyloseq"))
#' unlink(paste0(tempdir(), "/phyloseq"), recursive = TRUE)
read_pq <- function(path = NULL,
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' See https://www.nature.com/articles/s41467-017-01312-x for more information
#'  on the method.
#'
#' @inheritParams clean_pq
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param id (default: 0.84) id for --usearch_global.
#' @param vsearchpath (default: vsearch) path to vsearch.
#' @param verbose (logical) If true, print some additional messages.
#' @param clean_pq (logical) If true, empty samples and empty ASV are discarded
#'   before clustering.
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#' @param ... Additional arguments passed on to function [lulu()]
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
#' @seealso [mumu_pq()]
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
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
                    keep_temporary_files = FALSE,
                    ...) {
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
    lulu(data.frame(t(physeq@otu_table), ...), match_list)

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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
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
#' @param vsearchpath (default: vsearch) path to vsearch.
#' @param mumupath path to mumu. See [mumu](https://github.com/frederic-mahe/mumu)
#'   for installation instruction
#' @param lulu_exact (logical) If true, use the exact same algorithm as LULU 
#'  corresponding to the --legacy option of mumu. Need mumu version >= v1.1.0
#' @param verbose (logical) If true, print some additional messages.
#' @param clean_pq (logical) If true, empty samples and empty ASV are discarded
#'   before clustering.
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#' @param extra_mumu_args (character, default: NULL) Additional arguments passed
#'  on to mumu command line. See `man mumu` into bash for details. Major args are
#' `--minimum_match`, `--minimum_ratio_type`, `--minimum_ratio`, 
#' `--minimum_relative_cooccurence` and `--threads` 
#' @return a list of for object
#' - "new_physeq": The new phyloseq object (class physeq)
#' - "mumu_results": The log file of the mumu software. Run `man mumu` into
#'   bash to obtain details about columns' signification.
#'
#' @export
#' @seealso [lulu_pq()]
#' @examplesIf MiscMetabar::is_mumu_installed()
#' \dontrun{
#' ntaxa(data_fungi_sp_known)
#' ntaxa(mumu_pq(data_fungi_sp_known)$new_physeq)
#' ntaxa(mumu_pq(data_fungi_sp_known, extra_mumu_args="--minimum_match 90")$new_physeq)
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
mumu_pq <- function(physeq,
                    nproc = 1,
                    id = 0.84,
                    vsearchpath = "vsearch",
                    mumupath = "mumu",
                    lulu_exact = FALSE,
                    verbose = FALSE,
                    clean_pq = TRUE,
                    keep_temporary_files = FALSE,
                    extra_mumu_args = NULL) {
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
    data.frame(unclass(taxa_as_rows(physeq)@otu_table))
  otu_tab <- cbind("Taxa" = rownames(otu_tab), otu_tab)
  write.table(
    otu_tab,
    "otu_table.csv",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  message("Mumu algorithm")

    mumu_cmd <-
    paste0( 
      " --otu_table otu_table.csv ",
      " --match_list match_list.txt ",
      " --log log.txt ",
      " --new_otu_table new_OTU.tablemumu"
    )
  
  if(!is.null(extra_mumu_args)){
    mumu_cmd <- paste0(mumu_cmd, " ", extra_mumu_args)
    message("mumu is runned with option(s)", extra_mumu_args)
  }
  if (lulu_exact) {
    mumu_cmd <- paste0(mumu_cmd, " --legacy ")
    message("Using LULU exact mode (--legacy option in mumu)")
  }
  system2(
    mumupath, args=mumu_cmd
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

  return(list("new_physeq" = new_physeq, "mumu_results" = result_mumu))
}
################################################################################




################################################################################
#' Verify the validity of a phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Mostly for internal use in MiscMetabar functions.
#'
#' @inheritParams clean_pq
#' @param verbose (logical, default FALSE) If TRUE, prompt some warnings.
#' @param min_nb_seq_sample (numeric) Only used if verbose = TRUE.
#'   Minimum number of sequences per samples to not show warning.
#' @param min_nb_seq_taxa (numeric) Only used if verbose = TRUE.
#'   Minimum number of sequences per taxa to not show warning.
#' @return Nothing if the phyloseq object is valid. An error in the other case.
#'  Warnings if verbose = TRUE
#' @export
#' @author Adrien Taudière
verify_pq <- function(physeq,
                      verbose = FALSE,
                      min_nb_seq_sample = 500,
                      min_nb_seq_taxa = 1) {
  if (!methods::validObject(physeq) ||
    !inherits(physeq, "phyloseq")) {
    stop("The physeq argument is not a valid phyloseq object.")
  }
  if (verbose) {
    if (min(sample_sums(physeq)) < min_nb_seq_sample) {
      warning(
        paste0(
          "At least one of your sample contains less than ",
          min_nb_seq_sample,
          " sequences."
        )
      )
    }
    if (min(sample_sums(physeq)) < min_nb_seq_sample) {
      warning(
        paste0(
          "At least one of your taxa is represent by less than ",
          min_nb_seq_taxa,
          " sequences."
        )
      )
    }
    if (sum(is.na(physeq@sam_data)) > 0) {
      warning("At least one of your samples metadata columns contains NA.")
    }
    if (sum(grepl("^[0-9]", sample_names(physeq)) > 0)) {
      message(
        "At least one sample name start with a number.
      It may introduce bug in some function such
      as psmelt. You may replace sample_names using
      for example :
      sample_names(physeq) <- paste('samp', sample_names(physeq))"
      )
    }
  }
}
################################################################################


################################################################################
#' Subset samples using a conditional boolean vector.
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
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
#' @author Adrien Taudière
#' @examples
#'
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
  if (is.null(physeq@sam_data)) {
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' The main objective of this function is to complete the
#' [phyloseq::subset_taxa()] function by propose a more easy way of
#' subset_taxa using a named boolean vector. Names must match taxa_names.
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
#'
#' subset_taxa_pq(data_fungi, data_fungi@tax_table[, "Phylum"] == "Ascomycota")
#' subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 100)
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
#' @author Adrien Taudière
subset_taxa_pq <- function(physeq,
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
      paste(names(condition)[!names(condition) %in% taxa_names(physeq)], collapse = "/")
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
  new_MA <- old_MA[cond, , drop = FALSE]

  new_otu_table <- otu_table(new_MA, taxa_are_rows = TRUE)

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
#' Filter taxa of a phyloseq object based on the minimum number of
#'   sequences/samples
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Basically a wraper of [subset_taxa_pq()].
#'
#' @inheritParams clean_pq
#' @param min_nb_seq (int default NULL) minimum number of sequences by taxa.
#' @param min_occurence (int default NULL) minimum number of sample by taxa.
#' @param combination Either "AND" (default) or "OR". If set to "AND" and both
#'   min_nb_seq and min_occurence are not NULL, the taxa must match the two
#'   condition to passe the filter. If set to "OR", taxa matching only one
#'   condition are kept.
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples and empty taxa (ASV, OTU) are discarded
#'   after filtering.
#'
#' @return a new phyloseq object
#' @export
#' @author Adrien Taudière
#' @examples
#' filt_taxa_pq(data_fungi, min_nb_seq = 20)
#' filt_taxa_pq(data_fungi, min_occurence = 2)
#' filt_taxa_pq(data_fungi,
#'   min_occurence = 2,
#'   min_nb_seq = 10, clean_pq = FALSE
#' )
#' filt_taxa_pq(data_fungi,
#'   min_occurence = 2,
#'   min_nb_seq = 10,
#'   combination = "OR"
#' )
filt_taxa_pq <- function(physeq,
                         min_nb_seq = NULL,
                         min_occurence = NULL,
                         combination = "AND",
                         clean_pq = TRUE) {
  new_physeq <- physeq
  if (is.null(min_nb_seq) && is.null(min_occurence)) {
    stop("You must inform either min_nb_seq or min_occurence or both params!")
  }
  if (!is.null(min_nb_seq) && is.null(min_occurence)) {
    new_physeq <-
      subset_taxa_pq(physeq, taxa_sums(physeq) >= min_nb_seq)
  }

  if (is.null(min_nb_seq) && !is.null(min_occurence)) {
    new_physeq <-
      subset_taxa_pq(
        new_physeq,
        taxa_sums(as_binary_otu_table(new_physeq)) >= min_occurence
      )
  }

  if (!is.null(min_nb_seq) && !is.null(min_occurence)) {
    if (combination == "AND") {
      new_physeq <-
        subset_taxa_pq(
          new_physeq,
          taxa_sums(as_binary_otu_table(new_physeq)) >= 1 &
            taxa_sums((new_physeq)) >= 100
        )
    } else if (combination == "OR") {
      new_physeq <-
        subset_taxa_pq(
          new_physeq,
          taxa_sums(as_binary_otu_table(new_physeq)) >= 1 |
            taxa_sums((new_physeq)) >= 100
        )
    }
  }

  if (clean_pq) {
    new_physeq <- clean_pq(new_physeq)
  }
  return(new_physeq)
}


################################################################################
#' Select one sample from a physeq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Mostly for internal used, for example in function [track_wkflow_samples()].
#'
#' @inheritParams clean_pq
#' @param sam_name (required) The sample name to select
#' @param silent (logical) If true, no message are printing.
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with one sample
#'
#' @export
#'
#' @author Adrien Taudière
#'
#' @examples
#'
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' One of main use of this function is to add taxonomic assignment from
#' a new database.
#'
#' @inheritParams clean_pq
#' @param ref_fasta (required) A link to a database.
#'   passed on to `dada2::assignTaxonomy`.
#' @param method (required, default "dada2") :
#'
#' - "dada2": [dada2::assignTaxonomy()]
#'
#' - "dada2_2step": [assign_dada2()]
#'
#' - "sintax": see [assign_sintax()]
#'
#' - "lca": see [assign_vsearch_lca()]
#'
#' - "idtaxa": see [assign_idtaxa()]
#'
#' - "blastn": see [assign_blastn()]
#'
#' @param suffix (character) The suffix to name the new columns.
#'   If set to NULL (the default), the basename of the file reFasta
#'   is used with the name of the method. Set suffix to "" in order
#'   to remove any suffix.
#' @param trainingSet see [assign_idtaxa()]. Only used if method = "idtaxa".
#'   Note that if trainingSet is not NULL, the ref_fasta is overwrite by the
#'   trainingSet parameter. To customize learning parameters of the idtaxa
#'   algorithm you must use trainingSet computed by the function [learn_idtaxa()].
#' @param min_bootstrap  (Float \[0:1\])
#'
#'   Minimum bootstrap value to inform taxonomy. For each bootstrap
#'   below the min_bootstrap value, the taxonomy information is set to NA.
#'
#'   Correspond to parameters :
#'
#'  - dada2 & dada2_2step: `minBoot`, default value = 0.5
#'
#'  - sintax: `min_bootstrap`, default value = 0.5
#'
#'  - lca: `id`, default value = 0.5. Note in that case, the bootstrap value is different.
#'    See the id parameter in [assign_vsearch_lca()]
#'
#'  - idtaxa: `threshold`, default value = 0.6
#'
#'  - blastn: This method do not take different bootstrap value. You may use method="vote" with
#'    different `vote_algorithm` as well as different filters parameters
#'    (min_id, min_bit_score, min_cover and min_e_value)
#'
#' @param ... Additional arguments passed on to the taxonomic assignation method.
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a larger slot tax_table"
#' @seealso [dada2::assignTaxonomy()], [assign_sintax()], [assign_vsearch_lca()], [assign_sintax()], [assign_blastn()], [assign_dada2()]
#' @export
#' @examples
#' \dontrun{
#' ref_fasta <- system.file("extdata",
#'   "mini_UNITE_fungi.fasta.gz",
#'   package = "MiscMetabar", mustWork = TRUE
#' )
#' add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "dada2")
#' add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "dada2_2steps")
#' add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "lca")
#' add_new_taxonomy_pq(data_fungi_mini, ref_fasta, method = "idtaxa")
#'
#' # blastn doesn't work with fasta.gz format
#' ref_fasta <- system.file("extdata",
#'   "100_sp_UNITE_sh_general_release_dynamic_sintax.fasta",
#'   package = "MiscMetabar", mustWork = TRUE
#' )
#'
#' dp <- add_new_taxonomy_pq(data_fungi_mini, ref_fasta,
#'   method = "blastn", min_id = 80, min_cover = 50, min_bit_score = 20,
#'   min_e_value = 1e-20
#' )
#' dp_tophit <- add_new_taxonomy_pq(data_fungi_mini, ref_fasta,
#'   method = "blastn", min_id = 80, min_cover = 50, min_bit_score = 20,
#'   min_e_value = 1e-20, method_algo = "top_hit"
#' )
#' }
#' @author Adrien Taudière
#'
add_new_taxonomy_pq <- function(
    physeq,
    ref_fasta,
    suffix = NULL,
    method = c("dada2", "sintax", "lca", "idtaxa", "blastn", "dada2_2steps"),
    trainingSet = NULL,
    min_bootstrap = NULL,
    ...) {
  method <- match.arg(method)

  if (is.null(min_bootstrap)) {
    min_bootstrap <- ifelse(method == "idtaxa", 0.6, 0.5)
  }

  if (is.null(suffix)) {
    suffix <- paste0("_", basename(ref_fasta), "_", method)
  }
  if (method == "dada2") {
    list_args <- list(seqs = physeq@refseq, refFasta = ref_fasta, minBoot = 100 * min_bootstrap, ...)
    list_args <- list_args[names(list_args) %in% c(
      "seqs", "refFasta", "minBoot", "tryRC", "outputBootstraps",
      "taxLevels", "multithread", "verbose"
    )]
    tax_tab <- do.call(dada2::assignTaxonomy, args = list_args, envir = parent.frame())
    colnames(tax_tab) <-
      make.unique(paste0(colnames(tax_tab), suffix))
    new_tax_tab <- tax_table(cbind(physeq@tax_table, tax_tab))
    new_physeq <- physeq
    tax_table(new_physeq) <- new_tax_tab
  } else if (method == "dada2_2steps") {
    list_args <- list(physeq = physeq, ref_fasta = ref_fasta, min_bootstrap = min_bootstrap, suffix = suffix, ...)
    list_args <- list_args[names(list_args) %in% names(formals("assign_dada2", envir = parent.frame()))]
    new_physeq <- do.call(assign_dada2, args = list_args, envir = parent.frame())
  } else if (method == "sintax") {
    list_args <- list(physeq = physeq, ref_fasta = ref_fasta, suffix = suffix, behavior = "add_to_phyloseq", min_bootstrap = min_bootstrap, ...)
    list_args <- list_args[names(list_args) %in% names(formals("assign_sintax", envir = parent.frame()))]
    new_physeq <- do.call(assign_sintax, args = list_args, envir = parent.frame())
  } else if (method == "lca") {
    list_args <- list(physeq = physeq, ref_fasta = ref_fasta, suffix = suffix, behavior = "add_to_phyloseq", ...)
    list_args <- list_args[names(list_args) %in% names(formals("assign_vsearch_lca", envir = parent.frame()))]
    new_physeq <- do.call(assign_vsearch_lca, args = list_args, envir = parent.frame())
  } else if (method == "idtaxa") {
    if (is.null(trainingSet)) {
      list_args <- list(physeq = physeq, fasta_for_training = ref_fasta, suffix = suffix, behavior = "add_to_phyloseq", threshold = 100 * min_bootstrap, ...)
      list_args <- list_args[names(list_args) %in% names(formals("assign_idtaxa", envir = parent.frame()))]
      new_physeq <- do.call(assign_idtaxa, args = list_args, envir = parent.frame())
    } else {
      list_args <- list(physeq = physeq, trainingSet = trainingSet, suffix = suffix, behavior = "add_to_phyloseq", threshold = 100 * min_bootstrap, ...)
      list_args <- list_args[names(list_args) %in% names(formals("assign_idtaxa", envir = parent.frame()))]
      new_physeq <- do.call(assign_idtaxa, args = list_args, envir = parent.frame())
    }
  } else if (method == "blastn") {
    list_args <- list(physeq = physeq, ref_fasta = ref_fasta, suffix = suffix, behavior = "add_to_phyloseq", ...)
    list_args <- list_args[names(list_args) %in% names(formals("assign_blastn", envir = parent.frame()))]
    new_physeq <- do.call(assign_blastn, args = list_args, envir = parent.frame())
  }

  return(new_physeq)
}
################################################################################


################################################################################
#' Summarize information from sample data in a table
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A wrapper for the [gtsummary::tbl_summary()] function in the case of `physeq`
#'   object.
#'
#' @inheritParams clean_pq
#' @param remove_col_unique_value (logical, default TRUE) Do we remove
#'  informative columns (categorical column with one value per samples),
#'   e.g. samples names ?
#' @param ... Additional arguments passed on to [gtsummary::tbl_summary()].
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with a larger slot tax_table
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' if (requireNamespace("gtsummary")) {
#'   tbl_sum_samdata(data_fungi) %>%
#'     gtsummary::as_kable()
#'
#'   summary_samdata <- tbl_sum_samdata(data_fungi,
#'     include = c("Time", "Height"),
#'     type = list(Time ~ "continuous2", Height ~ "categorical"),
#'     statistic = list(Time ~ c("{median} ({p25}, {p75})", "{min}, {max}"))
#'   )
#' }
#' \donttest{
#' data(enterotype)
#' if (requireNamespace("gtsummary")) {
#'   summary_samdata <- tbl_sum_samdata(enterotype)
#'   summary_samdata <- tbl_sum_samdata(enterotype, include = !contains("SampleId"))
#' }
#' }
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

################################################################################
#' Summarize a tax_table (taxonomic slot of phyloseq object) using gtsummary
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Mainly a wrapper for the [gtsummary::tbl_summary()] function in the case
#' of `physeq` object.
#'
#' @inheritParams clean_pq
#' @param taxonomic_ranks A list of taxonomic ranks we want to summarized.
#' @param ... Additional arguments passed on to [gtsummary::tbl_summary()]
#'
#' @return A table of class c('tbl_summary', 'gtsummary')
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' tbl_sum_taxtable(data_fungi_mini)
#' data_fungi_mini |>
#'   filt_taxa_pq(min_occurence = 2) |>
#'   tbl_sum_taxtable(taxonomic_rank = c("Species", "Genus"))
tbl_sum_taxtable <- function(physeq, taxonomic_ranks = NULL, ...) {
  taxatab <- as.data.frame(unclass(physeq@tax_table))

  if (!is.null(taxonomic_ranks)) {
    taxatab <- taxatab |>
      select(one_of(taxonomic_ranks))
  }
  tbl_sum <- gtsummary::tbl_summary(taxatab, ...)
  return(tbl_sum)
}
################################################################################

################################################################################
#' Add information about Guild for FUNGI the FUNGuild databse
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Please cite Nguyen et al. 2016 (\doi{doi:10.1016/j.funeco.2015.06.006})
#'
#' @inheritParams clean_pq
#' @param taxLevels Name of the 7 columns in tax_table required by funguild
#' @param db_url a length 1 character string giving the URL to retrieve the database
#'     from
#' @return A new object of class `physeq` with Guild information added to
#'   `tax_table` slot
#' @export
#' @author Adrien Taudière
#' @examples
#' \dontrun{
#' # to avoid bug in CRAN when internet is not available
#' if (requireNamespace("httr")) {
#'   d_fung_mini <- add_funguild_info(data_fungi_mini,
#'     taxLevels = c(
#'       "Domain",
#'       "Phylum",
#'       "Class",
#'       "Order",
#'       "Family",
#'       "Genus",
#'       "Species"
#'     )
#'   )
#'   sort(table(d_fung_mini@tax_table[, "guild"]), decreasing = TRUE)
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `FUNGuildR` package and the associate
#'   publication (\doi{doi:10.1016/j.funeco.2015.06.006}) if you
#'   use this function.
#' @seealso [plot_guild_pq()]

add_funguild_info <- function(physeq,
                              taxLevels = c(
                                "Kingdom",
                                "Phylum",
                                "Class",
                                "Order",
                                "Family",
                                "Genus",
                                "Species"
                              ),
                              db_url = "http://www.stbates.org/funguild_db_2.php") {
  if (httr::http_error(db_url)) {
    message("error with db_url: ", db_url)
    return(NULL)
  }

  tax_tab <- physeq@tax_table
  FUNGuild_assign <-
    funguild_assign(data.frame(
      "Taxonomy" =
        apply(tax_tab[, taxLevels], 1, paste, collapse = ";")
    ), db_url = db_url)
  if (is.null(FUNGuild_assign)) {
    message("No http access to the funguild database. No information were added.")
    return(physeq)
  }
  tax_tab <-
    as.matrix(cbind(tax_tab, FUNGuild_assign))
  physeq@tax_table <- tax_table(tax_tab)
  return(physeq)
}




################################################################################
#' Plot information about Guild from tax_table slot previously
#' created with [add_funguild_info()]
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Graphical function.
#'
#' @inheritParams clean_pq
#' @param levels_order (Default NULL) A character vector to
#'   reorder the levels of guild. See examples.
#' @param clean_pq (logical, default TRUE): Does the phyloseq
#'   object is cleaned using the [clean_pq()] function?
#' @param ... Other params for be passed on to
#'   [clean_pq()] function
#' @return A ggplot2 object
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' \dontrun{
#' # to avoid bug in CRAN when internet is not available
#' if (requireNamespace("httr")) {
#'   d_fung_mini <- add_funguild_info(data_fungi_mini,
#'     taxLevels = c(
#'       "Domain",
#'       "Phylum",
#'       "Class",
#'       "Order",
#'       "Family",
#'       "Genus",
#'       "Species"
#'     )
#'   )
#'   sort(table(d_fung_mini@tax_table[, "guild"]), decreasing = TRUE)
#'
#'   p <- plot_guild_pq(d_fung_mini)
#'   if (requireNamespace("patchwork")) {
#'     (plot_guild_pq(subset_samples(d_fung_mini, Height == "Low"),
#'       levels_order = p$data$Guild[order(p$data$nb_seq)]
#'     ) + theme(legend.position = "none")) +
#'       (plot_guild_pq(subset_samples(d_fung_mini, Height == "High"),
#'         levels_order = p$data$Guild[order(p$data$nb_seq)]
#'       ) + ylab("") + theme(axis.text.y = element_blank()))
#'   }
#' }
#' }
#' @seealso [add_funguild_info()]

plot_guild_pq <-
  function(physeq,
           levels_order = NULL,
           clean_pq = TRUE,
           ...) {
    if (clean_pq) {
      physeq <- clean_pq(physeq, ...)
    }
    guilds <-
      data.frame(sort(table(strsplit(
        paste(physeq@tax_table[, "guild"]
        [physeq@tax_table[, "confidenceRanking"] %in%
            c("Highly Probable", "Probable")], collapse = "-"),
        split = "-"
      ))))

    guilds$Var1 <- as.vector(guilds$Var1)
    guilds <- guilds[guilds$Var1 != "NA", ]
    guilds <- guilds[guilds$Var1 != "NULL", ]
    guilds <- guilds[guilds$Var1 != "", ]

    # Number of sequences per guild
    nb_seq_by_guild <- vector("integer", length(guilds$Var1))
    for (i in seq_along(guilds$Var1)) {
      nb_seq_by_guild[i] <-
        sum(taxa_sums(physeq@otu_table)[grepl(guilds$Var1[i], physeq@tax_table[, "guild"])])
    }
    names(nb_seq_by_guild) <- guilds$Var1
    guilds$seq <- nb_seq_by_guild

    names(guilds) <- c("Guild", "nb_taxa", "nb_seq")
    guilds$nb_seq <- as.numeric(guilds$nb_seq)
    guilds$nb_taxa <- as.numeric(guilds$nb_taxa)

    guilds$Guild <- factor(as.vector(guilds$Guild), levels = guilds$Guild[order(guilds$nb_seq)])


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
        "nb_taxa" = ntaxa(physeq),
        "nb_seq" = sum(physeq@otu_table),
        "colors" = "ALL",
        stringsAsFactors = FALSE
      )
    )
    guilds <- guilds[order(guilds$nb_seq), ]
    if (!is.null(levels_order)) {
      guilds$Guild <- factor(guilds$Guild, levels = levels_order)
    }

    ggplot(guilds, aes(
      y = Guild,
      x = log10(nb_seq),
      fill = colors
    )) +
      geom_bar(stat = "identity") +
      annotation_logticks(sides = "b", alpha = 0.5) +
      ylab("GUILD by FUNGuild") +
      scale_fill_manual("Guild",
        values = c("gray", "Olivedrab", "cyan4", "tomato3", "lightpink4")
      ) +
      geom_text(aes(label = nb_taxa, x = log10(nb_seq) + 0.2), family = "serif") +
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' This function build tree phylogenetic tree and if nb_bootstrap is
#' set, it build also the 3 corresponding bootstrapped tree.
#'
#' Default parameters are based on \doi{doi:10.12688/f1000research.8986.2}
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
#' @param ... Other params for be passed on to
#'   [phangorn::optim.pml()] function
#'
#' @return A list of phylogenetic tree
#' @export
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `phangorn` package if you
#'   use this function.
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' \donttest{
#' if (requireNamespace("phangorn")) {
#'   set.seed(22)
#'   df <- subset_taxa_pq(data_fungi_mini, taxa_sums(data_fungi_mini) > 9000)
#'   df_tree <- build_phytree_pq(df, nb_bootstrap = 2)
#'   plot(df_tree$UPGMA)
#'   phangorn::plotBS(df_tree$UPGMA, df_tree$UPGMA_bs, main = "UPGMA")
#'   plot(df_tree$NJ, "unrooted")
#'   plot(df_tree$ML)
#'
#'   phangorn::plotBS(df_tree$ML$tree, df_tree$ML_bs, p = 20, frame = "circle")
#'   phangorn::plotBS(
#'     df_tree$ML$tree,
#'     df_tree$ML_bs,
#'     p = 20,
#'     frame = "circle",
#'     method = "TBE"
#'   )
#'   plot(phangorn::consensusNet(df_tree$ML_bs))
#'   plot(phangorn::consensusNet(df_tree$NJ_bs))
#'   ps_tree <- merge_phyloseq(df, df_tree$ML$tree)
#' }
#' }
build_phytree_pq <- function(physeq,
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
      phangorn::bootstrap.phyDat(phang.align, function(x) {
        phangorn::upgma(phangorn::dist.ml(x))
      }, bs = nb_bootstrap)
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
    treeNJ_bs <- phangorn::bootstrap.phyDat(phang.align, function(x) {
      phangorn::NJ(phangorn::dist.ml(x))
    }, bs = nb_bootstrap)
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' The aim of this function is to provide a warnings if samples depth significantly
#' vary among the modalities of a factor present in the `sam_data` slot.
#'
#' This function apply a Kruskal-Wallis rank sum test to the number of sequences
#' per samples in function of the factor `fact`.
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}.
#' @param boxplot (logical) Do you want to plot boxplot?
#'
#' @return The result of a Kruskal-Wallis rank sum test
#' @export
#' @author Adrien Taudière
#' @importFrom stats kruskal.test
#' @examples
#'
#' are_modality_even_depth(data_fungi_mini, "Time")$p.value
#' are_modality_even_depth(rarefy_even_depth(data_fungi_mini), "Time")$p.value
#' are_modality_even_depth(data_fungi_mini, "Height", boxplot = TRUE)
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Note that the taxa order in a physeq object with a tree is locked by
#' the order of leaf in the phylogenetic tree.
#'
#' @inheritParams clean_pq
#' @param names_ordered (required) Names of the taxa (must be the same
#'   as taxa in `taxa_names(physeq)`) in a given order
#' @param remove_phy_tree (logical, default FALSE) If TRUE, the phylogenetic
#'   tree is removed. It is
#' @return A phyloseq object
#' @export
#' @author Adrien Taudière
#' @examples
#'
#' data_fungi_ordered_by_genus <- reorder_taxa_pq(
#'   data_fungi,
#'   taxa_names(data_fungi)[order(as.vector(data_fungi@tax_table[, "Genus"]))]
#' )
#'
#' data_fungi_mini_asc_ordered_by_abundance <- reorder_taxa_pq(
#'   data_fungi_mini,
#'   taxa_names(data_fungi_mini)[order(taxa_sums(data_fungi_mini))]
#' )
reorder_taxa_pq <- function(physeq, names_ordered, remove_phy_tree = FALSE) {
  new_physeq <- physeq

  if (!is.null(phy_tree(new_physeq, FALSE))) {
    if (remove_phy_tree) {
      message("Removing phylogenetic tree!")
      new_physeq@phy_tree <- NULL
    } else {
      stop(
        "The taxa order in a physeq object with a tree is locked by
      the order of leaf in the phylogenetic tree. You could use args
      remove_phy_tree = TRUE."
      )
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

################################################################################
#' @title Add information to sample_data slot of a phyloseq-class object
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Warning: The value nb_seq and nb_otu may be outdated if you transform your
#' phyloseq object, e.g. using the [subset_taxa_pq()] function
#'
#' @inheritParams clean_pq
#' @param df_info : A dataframe with rownames matching for sample names of the
#'   phyloseq object
#' @param add_nb_seq (Logical, default TRUE) Does we add a column nb_seq
#'   collecting the number of sequences per sample?
#' @param add_nb_otu (Logical, default TRUE) Does we add a column nb_otu
#'   collecting the number of OTUs per sample?
#'
#' @return A phyloseq object with an updated sam_data slot
#' @export
#'
#' @examples
#'
#' data_fungi <- add_info_to_sam_data(data_fungi)
#' boxplot(data_fungi@sam_data$nb_otu ~ data_fungi@sam_data$Time)
#'
#' new_df <- data.frame(
#'   variable_1 = runif(n = nsamples(data_fungi), min = 1, max = 20),
#'   variable_2 = runif(n = nsamples(data_fungi), min = 1, max = 2)
#' )
#' rownames(new_df) <- sample_names(data_fungi)
#' data_fungi <- add_info_to_sam_data(data_fungi, new_df)
#' plot(data_fungi@sam_data$nb_otu ~ data_fungi@sam_data$variable_1)
#' @author Adrien Taudière

add_info_to_sam_data <- function(physeq,
                                 df_info = NULL,
                                 add_nb_seq = TRUE,
                                 add_nb_otu = TRUE) {
  if (add_nb_seq) {
    physeq@sam_data$nb_seq <- sample_sums(physeq)
  }
  if (add_nb_otu) {
    physeq@sam_data$nb_otu <- sample_sums(as_binary_otu_table(physeq))
  }
  if (!is.null(df_info)) {
    if (sum(sample_names(physeq) %in% rownames(df_info)) == 0) {
      stop("Rownames of df_info must match the sample names of physeq.")
    }
    df_info <-
      df_info[match(sample_names(physeq), rownames(df_info)), ]
    physeq@sam_data <- sample_data(cbind(as.data.frame(physeq@sam_data), df_info))
  }
  return(physeq)
}
################################################################################


###############################################################################
#' Return a DNAStringSet object from either a character vector of DNA sequences
#'   or the `refseq` slot of a phyloseq-class object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#'   Internally used in [vsearch_clustering()], [swarm_clustering()] and
#'   [postcluster_pq()].
#'
#' @inheritParams clean_pq
#' @param dna_seq You may directly use a character vector of DNA sequences
#'   in place of physeq args. When physeq is set, dna sequences take the value
#'   of `physeq@refseq`
#'
#' @return An object of class DNAStringSet (see the [Biostrings::DNAStringSet()]
#'   function)
#' @export
#'
#' @examples
#'
#' dna <- physeq_or_string_to_dna(data_fungi)
#' dna
#'
#' sequences_ex <- c(
#'   "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTAATAACGAATTCATTGAATCA",
#'   "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAAGGCCCACTT",
#'   "TACCTATGTTGCCTTGGCGGCTAAACCTACCCGGGATTTGATGGGGCGAATTACCTGGTAGAGGTG",
#'   "TACCTATGTTGCCTTGGCGGCTAAACCTACC",
#'   "CGGGATTTGATGGCGAATTACCTGGTATTTTAGCCCACTTACCCGGTACCATGAGGTG",
#'   "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACCTGG",
#'   "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
#'   "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG",
#'   "GCGGCTAAACCTACCCGGGATTTGATGGCGAATTACAAAG"
#' )
#' dna2 <- physeq_or_string_to_dna(dna_seq = sequences_ex)
#' dna2
#'
#' @seealso [Biostrings::DNAStringSet()]
#' @author Adrien Taudière
physeq_or_string_to_dna <- function(physeq = NULL, dna_seq = NULL) {
  if (inherits(physeq, "phyloseq")) {
    verify_pq(physeq)
    if (is.null(physeq@refseq)) {
      stop("The phyloseq object do not contain a @refseq slot")
    }
    dna <- Biostrings::DNAStringSet(physeq@refseq)
    if (!is.null(dna_seq)) {
      stop("You must use either physeq or dna_seq args but not both")
    }
  } else if (inherits(dna_seq, "character")) {
    dna <- Biostrings::DNAStringSet(dna_seq)
  } else if (inherits(dna_seq, "DNAStringSet")) {
    return(dna)
  } else {
    stop(
      "You must set the args physeq (object of class phyloseq) or
    dna_seq (character vector)."
    )
  }
  return(dna)
}
###############################################################################



################################################################################
#' Remove primers using [cutadapt](https://github.com/marcelm/cutadapt/)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' You need to install [Cutadapt](https://cutadapt.readthedocs.io/).
#'   See also https://github.com/VascoElbrecht/JAMP/blob/master/JAMP/R/Cutadapt.R for another call to cutadapt
#'   from R
#'
#' @param path_to_fastq (Required) A path to a folder with fastq files. See
#'   [list_fastq_files()] for help.
#' @inheritParams list_fastq_files
#' @param primer_fw (Required, String) The forward primer DNA sequence.
#' @param primer_rev (String)  The reverse primer DNA sequence.
#' @param folder_output The path to a folder for output files
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param cmd_is_run (logical, default TRUE) Do the cutadapt command is run.
#'   If set to FALSE, the only effect of the function is to return a list of
#'   command to manually run in a terminal.
#' @param return_file_path (logical, default FALSE) If true, the function
#'   return the path of the output folder (param `folder_output`). Useful
#'   in targets workflow
#' @param args_before_cutadapt (String) A one line bash command to run before
#' to run cutadapt. For examples, "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv &&" allow to bypass the conda init which asks to restart the shell
#'
#' @return a list of command or if `return_file_path` is TRUE, the path to
#'   the output folder
#' @export
#' @author Adrien Taudière
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' \dontrun{
#' cutadapt_remove_primers(system.file("extdata", package = "MiscMetabar"),
#'   "TTC",
#'   "GAA",
#'   folder_output = tempdir()
#' )
#'
#' cutadapt_remove_primers(
#'   system.file("extdata",
#'     package = "dada2"
#'   ),
#'   pattern_R1 = "F.fastq.gz",
#'   pattern_R2 = "R.fastq.gz",
#'   primer_fw = "TTC",
#'   primer_rev = "GAA",
#'   folder_output = tempdir()
#' )
#'
#' cutadapt_remove_primers(
#'   system.file("extdata",
#'     package = "dada2"
#'   ),
#'   pattern_R1 = "F.fastq.gz",
#'   primer_fw = "TTC",
#'   folder_output = tempdir(),
#'   cmd_is_run = FALSE
#' )
#'
#' unlink(tempdir(), recursive = TRUE)
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite cutadapt (\doi{doi:10.14806/ej.17.1.200}).


cutadapt_remove_primers <- function(path_to_fastq,
                                    primer_fw = NULL,
                                    primer_rev = NULL,
                                    folder_output = "wo_primers",
                                    nproc = 1,
                                    pattern = "fastq.gz",
                                    pattern_R1 = "_R1",
                                    pattern_R2 = "_R2",
                                    nb_files = Inf,
                                    cmd_is_run = TRUE,
                                    return_file_path = FALSE,
                                    args_before_cutadapt = "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv && ") {
  if (!dir.exists(folder_output)) {
    dir.create(folder_output)
  }


  if (is.null(primer_rev)) {
    lff <- list_fastq_files(
      path_to_fastq,
      paired_end = FALSE,
      pattern = pattern,
      pattern_R1 = pattern_R1,
      pattern_R2 = pattern_R2,
      nb_files = nb_files
    )
    cmd <- vector("list", length(lff$fnfs)) # pre-allocated for performance
    names(cmd) <- lff$fnfs
    for (f in lff$fnfs) {
      cmd[[f]] <-
        paste0(
          args_before_cutadapt,
          "cutadapt --cores=",
          nproc,
          " --json=",
          folder_output,
          "/",
          # More efficient: remove both extensions in one pass
          sub("\\.fastq(\\.gz)?$", "", basename(f)),
          ".cutadapt.json",
          " --discard-untrimmed -g '",
          primer_fw,
          "' -o ",
          folder_output,
          "/",
          basename(f),
          " ",
          f
        )
    }
  } else {
    lff <- list_fastq_files(
      path_to_fastq,
      paired_end = TRUE,
      pattern = pattern,
      pattern_R1 = pattern_R1,
      pattern_R2 = pattern_R2,
      nb_files = nb_files
    )
    cmd <- vector("list", length(lff$fnfs)) # pre-allocated for performance
    names(cmd) <- lff$fnfs

    primer_fw_RC <- dada2::rc(primer_fw)
    primer_rev_RC <- dada2::rc(primer_rev)
    for (f in lff$fnfs) {
      cmd[[f]] <-
        paste0(
          args_before_cutadapt,
          "cutadapt -n 2 --cores=",
          nproc,
          " --json=",
          folder_output,
          "/",
          # More efficient: remove both extensions in one pass
          sub("\\.fastq(\\.gz)?$", "", basename(f)),
          ".cutadapt.json",
          " --discard-untrimmed -g '",
          primer_fw,
          "' -G '",
          primer_rev,
          "' -a '",
          primer_rev_RC,
          "' -A '",
          primer_fw_RC,
          "' -o ",
          folder_output,
          "/",
          basename(f),
          " -p ",
          folder_output,
          "/",
          gsub(pattern_R1, pattern_R2, basename(f)),
          " ",
          f,
          " ",
          gsub(pattern_R1, pattern_R2, f)
        )
    }
  }
  if (cmd_is_run) {
    writeLines(unlist(cmd), paste0(tempdir(), "/script_cutadapt.sh"))
    system2("bash", paste0(tempdir(), "/script_cutadapt.sh"))
    message(paste0(
      "Output files are available in the folder ",
      normalizePath(folder_output)
    ))
    unlink(paste0(tempdir(), "/script_cutadapt.sh"))
  }
  if (return_file_path) {
    return(normalizePath(folder_output))
  } else {
    return(cmd)
  }
}
################################################################################

################################################################################
#' List the taxa founded only in one given level of a modality
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Given one modality name in sam_data and one level of the modality,
#' return the taxa strictly specific of this level.
#'
#' @inheritParams clean_pq
#' @param modality (required) The name of a column present in the `@sam_data` slot
#'   of the physeq object. Must be a character vector or a factor.
#' @param level (required) The level (must be present in modality) of interest
#' @param min_nb_seq_taxa (default 0 = no filter) The minimum number of sequences per taxa
#' @param min_nb_samples_taxa (default 0 = no filter) The minimum number of samples per taxa
#'
#' @return A vector of taxa names
#' @export
#'
#' @author Adrien Taudière
#' @examples
#' data_fungi_mini_woNA4height <- subset_samples(
#'   data_fungi_mini,
#'   !is.na(data_fungi_mini@sam_data$Height)
#' )
#' taxa_only_in_one_level(data_fungi_mini_woNA4height, "Height", "High")
#' #' # Taxa present only in low height samples
#' suppressMessages(suppressWarnings(
#'   taxa_only_in_one_level(data_fungi, "Height", "Low")
#' ))
#' # Number of taxa present only in sample of time equal to 15
#' suppressMessages(suppressWarnings(
#'   length(taxa_only_in_one_level(data_fungi, "Time", "15"))
#' ))
taxa_only_in_one_level <- function(physeq,
                                   modality,
                                   level,
                                   min_nb_seq_taxa = 0,
                                   min_nb_samples_taxa = 0) {
  if (min_nb_seq_taxa > 0) {
    physeq <-
      subset_taxa_pq(physeq, taxa_sums(physeq) >= min_nb_seq_taxa)
  }
  if (min_nb_samples_taxa > 0) {
    physeq <-
      subset_taxa_pq(
        physeq,
        taxa_sums(as_binary_otu_table(physeq)) >= min_nb_samples_taxa
      )
  }

  physeq_merged <- clean_pq(merge_samples2(physeq, modality))

  physeq_merged_only_one_level <-
    subset_taxa_pq(physeq_merged, taxa_sums(as_binary_otu_table(physeq_merged)) ==
      1)
  physeq_merged_only_level_given <-
    clean_pq(subset_samples_pq(
      physeq_merged_only_one_level,
      rownames(physeq_merged_only_one_level@sam_data) == level
    ))
  return(taxa_names(physeq_merged_only_level_given))
}
################################################################################


################################################################################
#' Normalize OTU table using samples depth
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   This function implement the method proposed by
#'   McKnight et al. 2018 (\doi{doi:10.5061/dryad.tn8qs35})
#'
#' @inheritParams clean_pq
#'
#' @param base_log (integer, default 2) the base for log-transformation. If
#'   set to NULL or NA, no log-transformation is compute after normalization.
#' @param constante a constante to multiply the otu_table values
#' @param digits (default = 2) integer indicating the number of decimal places
#'   to be used (see `?round` for more information)
#'
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object with otu_table count
#'   normalize and log transformed (if base_log is an integer)
#' @export
#' @author Adrien Taudière
#' @examples
#' taxa_sums(data_fungi_mini)
#' data_f_norm <- normalize_prop_pq(data_fungi_mini)
#' taxa_sums(data_f_norm)
#' ggplot(data.frame(
#'   "norm" = scale(taxa_sums(data_f_norm)),
#'   "raw" = scale(taxa_sums(data_fungi_mini)),
#'   "name_otu" = taxa_names(data_f_norm)
#' )) +
#'   geom_point(aes(x = raw, y = norm))
#'
#' data_f_norm <- normalize_prop_pq(data_fungi_mini, base_log = NULL)
normalize_prop_pq <- function(physeq,
                              base_log = 2,
                              constante = 10000,
                              digits = 4) {
  verify_pq(physeq)
  if (taxa_are_rows(physeq)) {
    new_otutab <- round((apply(physeq@otu_table, 2, function(x) {
      x / sum(x)
    })) * constante, digits = digits)
  } else {
    new_otutab <- round((apply(physeq@otu_table, 1, function(x) {
      x / sum(x)
    })) * constante, digits = digits)
  }

  if (!is.null(base_log) && !is.na(base_log)) {
    new_otutab <- round(log(new_otutab + 1, base = base_log), digits = digits)
  }

  new_physeq <- physeq
  new_physeq@otu_table <- otu_table(new_otutab, taxa_are_rows = taxa_are_rows(physeq))

  return(new_physeq)
}
################################################################################



################################################################################
#' Build a sample information tibble from physeq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Hill numbers are the number of equiprobable species giving the same diversity
#'   value as the observed distribution.
#'
#'   Note that contrary to [hill_pq()], this function does not take into
#'   account for difference in the number of sequences per samples/modalities.
#'   You may use rarefy_by_sample = TRUE if the mean number of sequences per
#'   samples differs among modalities.
#' @inheritParams clean_pq
#' @param hill_scales (a vector of integer) The list of q values to compute
#'   the hill number H^q. If Null, no hill number are computed. Default value
#'   compute the Hill number 0 (Species richness), the Hill number 1
#'   (exponential of Shannon Index) and the Hill number 2 (inverse of Simpson
#'   Index).
#' @param filter_zero (logical, default TRUE) Do we filter non present OTU from
#'   samples ? For the moment, this has no effect on the result because the dataframe
#'   is grouped by samples with abundance summed across OTU.
#' @param rarefy_by_sample (logical, default FALSE) If TRUE, rarefy
#'   samples using [phyloseq::rarefy_even_depth()] function.
#' @param rngseed (Optional). A single integer value passed to
#'   [phyloseq::rarefy_even_depth()], which is used to fix a seed for
#'   reproducibly random number generation (in this case, reproducibly
#'   random subsampling). If set to FALSE, then no fiddling with the RNG seed
#'   is performed, and it is up to the user to appropriately call set.seed
#'   beforehand to achieve reproducible results. Default is FALSE.
#' @param verbose (logical). If TRUE, print additional information.
#' @param taxa_ranks A vector of taxonomic ranks. For examples c("Family","Genus").
#'   If taxa ranks is not set (default value = NULL), taxonomic information are not
#'   present in the resulting tibble.
#' @author Adrien Taudière
#' @export
#' @return A tibble with a row for each sample. Columns provide information
#'   from `sam_data` slot as well as hill numbers, Abundance (nb of sequences),
#'   and Abundance_log10 (*log10(1+Abundance)*).
#' @examples
#' if (requireNamespace("ggstatsplot")) {
#'   psm_tib <- psmelt_samples_pq(data_fungi_mini, hill_scales = c(0, 2, 7))
#'   ggstatsplot::ggbetweenstats(psm_tib, Height, Hill_0)
#'   ggstatsplot::ggbetweenstats(psm_tib, Height, Hill_7)
#'
#'   psm_tib_tax <- psmelt_samples_pq(data_fungi_mini, taxa_ranks = c("Class", "Family"))
#'   ggplot(filter(psm_tib_tax, Abundance > 2000), aes(y = Family, x = Abundance, fill = Time)) +
#'     geom_bar(stat = "identity") +
#'     facet_wrap(~Height)
#' }
psmelt_samples_pq <-
  function(physeq,
           hill_scales = c(0, 1, 2),
           filter_zero = TRUE,
           rarefy_by_sample = FALSE,
           rngseed = FALSE,
           verbose = TRUE,
           taxa_ranks = NULL) {
    verify_pq(physeq)
    if (rarefy_by_sample) {
      if (as(rngseed, "logical")) {
        set.seed(rngseed)
        if (verbose) {
          message(
            "`set.seed(",
            rngseed,
            ")` was used to initialize repeatable random subsampling."
          )
          message("Please record this for your records so others can reproduce.")
          message("Try `set.seed(",
            rngseed,
            "); .Random.seed` for the full vector",
            sep = ""
          )
          message("...")
        }
      } else if (verbose) {
        message(
          "You set `rngseed` to FALSE. Make sure you've set & recorded\n",
          " the random seed of your session for reproducibility.\n",
          "See `?set.seed`\n"
        )
        message("...")
      }
      physeq <- rarefy_even_depth(physeq, rngseed = rngseed)
    }
    psm <- psmelt(physeq)
    if (filter_zero) {
      psm <- psm |> filter(Abundance > 0)
    }
    if (is.null(taxa_ranks)) {
      psm <- psm |>
        select(Sample, OTU, Abundance, colnames(physeq@sam_data))
      nb_distinct_samp <- psm |>
        group_by(Sample) |>
        select(-OTU, -Abundance) |>
        distinct() |>
        nrow()

      if (nsamples(physeq) != nb_distinct_samp) {
        stop(
          "The number of samples in physeq is different from the resulting
         number in psm tibble."
        )
      }
    } else {
      psm <- psm |>
        select(
          Sample,
          OTU,
          Abundance,
          colnames(physeq@sam_data),
          !!taxa_ranks
        )
    }

    if (is.null(taxa_ranks)) {
      psm_samp <- psm |>
        select(-OTU) |>
        group_by(Sample) |>
        summarise(
          Abundance = sum(Abundance),
          across(where(is.numeric) &
            !Abundance, ~ mean(.x, na.rm = TRUE)),
          across(where(is.character), ~ .x[1])
        )
    } else {
      psm_temp <- psm |>
        select(-OTU) |>
        group_by(Sample)

      for (i in seq_along(taxa_ranks)) {
        psm_temp <- psm_temp |>
          group_by(.data[[taxa_ranks[[i]]]], .add = TRUE)
      }

      psm_samp <- psm_temp |>
        summarise(
          Abundance = sum(Abundance),
          across(where(is.numeric) &
            !Abundance, ~ mean(.x, na.rm = TRUE)),
          across(where(is.character), ~ .x[1]),
          .groups = "drop"
        )
    }

    if (!is.null(hill_scales)) {
      physeq <- taxa_as_rows(physeq)
      df_hill <-
        vegan::renyi(t(physeq)@otu_table, scales = hill_scales, hill = TRUE)

      colnames(df_hill) <- paste0("Hill_", hill_scales)
      df_hill$Sample <- rownames(df_hill)

      psm_samp <- full_join(psm_samp, df_hill)
    }

    psm_samp <- psm_samp |> mutate(Abundance_log10 = log10(1 + Abundance))
    return(tibble(psm_samp))
  }
################################################################################

################################################################################
#' Force taxa to be in columns in the otu_table of a physeq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Mainly for internal use. It is a special case of clean_pq function.
#'
#' @inheritParams clean_pq
#' @author Adrien Taudière
#' @export
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object
taxa_as_columns <- function(physeq) {
  physeq <- clean_pq(
    physeq,
    clean_samples_names = FALSE,
    remove_empty_samples = FALSE,
    remove_empty_taxa = FALSE,
    force_taxa_as_columns = TRUE,
    silent = TRUE
  )
  return(physeq)
}
################################################################################

################################################################################
#' Force taxa to be in columns in the otu_table of a physeq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Mainly for internal use. It is a special case of clean_pq function.
#'
#' @inheritParams clean_pq
#' @author Adrien Taudière
#' @export
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object
taxa_as_rows <- function(physeq) {
  physeq <- clean_pq(
    physeq,
    clean_samples_names = FALSE,
    remove_empty_samples = FALSE,
    remove_empty_taxa = FALSE,
    force_taxa_as_rows = TRUE,
    silent = TRUE
  )
  return(physeq)
}
################################################################################

################################################################################
#' Rarefy (equalize) the number of samples per modality of a factor
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' This function randomly draw the same number of samples for each modality of factor.
#' It is usefull to dissentangle the effect of different number of samples per modality
#' on diversity. Internally used in [accu_plot_balanced_modality()].
#'
#' @inheritParams clean_pq
#' @param fact (required) The variable to rarefy. Must be present in
#'   the `sam_data` slot of the physeq object.
#' @param rngseed	(Optional). A single integer value passed to set.seed,
#'   which is used to fix a seed for reproducibly random number generation
#'   (in this case, reproducibly random subsampling). If set to FALSE, then no
#'   iddling with the RNG seed is performed, and it is up to the user to
#'   appropriately call
#' @param verbose (logical). If TRUE, print additional information.
#' @export
#' @author Adrien Taudière
#' @return A new \code{\link[phyloseq]{phyloseq-class}} object.
#' @seealso [accu_plot_balanced_modality()]
#' @examples
#' table(data_fungi_mini@sam_data$Height)
#' data_fungi_mini2 <- rarefy_sample_count_by_modality(data_fungi_mini, "Height")
#' table(data_fungi_mini2@sam_data$Height)
#' if (requireNamespace("patchwork")) {
#'   ggvenn_pq(data_fungi_mini, "Height") + ggvenn_pq(data_fungi_mini2, "Height")
#' }
rarefy_sample_count_by_modality <-
  function(physeq,
           fact,
           rngseed = FALSE,
           verbose = TRUE) {
    if (as(rngseed, "logical")) {
      set.seed(rngseed)
      if (verbose) {
        message(
          "`set.seed(",
          rngseed,
          ")` was used to initialize repeatable random subsampling."
        )
        message("Please record this for your records so others can reproduce.")
        message("Try `set.seed(",
          rngseed,
          "); .Random.seed` for the full vector",
          sep = ""
        )
        message("...")
      }
    } else if (verbose) {
      message(
        "You set `rngseed` to FALSE. Make sure you've set & recorded\n",
        " the random seed of your session for reproducibility.\n",
        "See `?set.seed`\n"
      )
      message("...")
    }
    mod <- as.factor(physeq@sam_data[[fact]])
    n_mod <- table(mod)
    samp_to_keep <- NULL
    for (modality in levels(mod)) {
      vec_samp_mod <- c(as.numeric(grep(modality, mod)))

      # To bypass the pb of vector of length 1
      # We build a vector of two equal values and we will take only one
      # It is cause by range base behavior:
      # 'If x has length 1, is numeric (in the sense of is.numeric) and x >= 1, sampling via sample takes place from 1:x.'
      if (length(vec_samp_mod) == 1) {
        vec_samp_mod <- c(vec_samp_mod, vec_samp_mod)
      }
      samp_to_keep <-
        c(
          samp_to_keep,
          sample(vec_samp_mod, size = min(n_mod), replace = FALSE)
        )
    }
    new_physeq <-
      subset_samples_pq(physeq, 1:nsamples(physeq) %in% samp_to_keep)

    if (length(table(new_physeq@sam_data[[fact]])) != length(table(mod))) {
      warning(
        paste0(
          "The number of final levels (sam_data of the output phyloseq
                    object) is not equal to the inital (sam_data of the input
                    phyloseq object) number of levels in the  factor: '",
          fact,
          "'"
        )
      )
    }

    return(new_physeq)
  }
################################################################################


################################################################################
#' A wrapper of \code{\link[DECIPHER]{IdTaxa}}
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' This function is basically a wrapper of functions [DECIPHER::IdTaxa()] and
#' [DECIPHER::LearnTaxa()], please cite the DECIPHER package if you use this
#' function. Note that if you want to specify parameters for the learning step
#' you must used the trainingSet param instead of the a fasta_for_training. The
#' training file can be obtain using the function [learn_idtaxa()].
#'
#' It requires:
#'
#' - either a physeq or seq2search object.
#'
#' - either a trainingSet or a fasta_for_training
#'
#' @inheritParams clean_pq
#' @param trainingSet An object of class Taxa and subclass
#'   Train compatible with the class of test.
#' @param seq2search A DNAStringSet object of sequences to search for. Replace
#'   the physeq object.
#' @param fasta_for_training A fasta file (can be gzip) to train the trainingSet
#'   using the function [learn_idtaxa()]. Only used if trainingSet is NULL.
#'
#'   The reference database must contain
#'   taxonomic information in the header of each sequence in the form of a string
#'   starting with ";tax=" and followed by a comma-separated list of up to nine
#'   taxonomic identifiers.
#'
#'   The only exception is if `unite=TRUE`. In that case the UNITE taxonomy is
#'   automatically formatted.
#' @param behavior Either "return_matrix" (default), or "add_to_phyloseq":
#'
#'  - "return_matrix" return a list of two objects. The first element is
#'    the taxonomic matrix and the second element is the raw results from
#'    DECIPHER::IdTaxa() function.
#'
#'  - "add_to_phyloseq" return a phyloseq object with amended slot `@taxtable`.
#'    Only available if using physeq input and not seq2search input.
#' @param threshold (Int, default 60) Numeric specifying the confidence at which
#'    to truncate the output taxonomic classifications.
#'   Lower values of threshold will classify deeper into the taxonomic tree at
#'   the expense of accuracy, and vise-versa for higher values of threshold. See
#'   [DECIPHER::IdTaxa()] man page.
#' @param column_names (vector of character) names for the column of the
#'   taxonomy
#' @param suffix (character) The suffix to name the new columns.
#'   Default to "_idtaxa".
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use
#' @param unite (logical, default FALSE). If set to TRUE, the fasta_for_training
#'   file is formatted from UNITE format to sintax one, needed in
#'   fasta_for_training. Only used if trainingSet is NULL.
#' @param verbose (logical). If TRUE, print additional information.
#' @param ... Additional arguments passed on to \code{\link[DECIPHER]{IdTaxa}}
#'
#' @seealso [assign_sintax()], [add_new_taxonomy_pq()], [assign_vsearch_lca()], [assign_blastn()]
#' @author Adrien Taudière
#' @return Either a new phyloseq object with additional information in
#'   the @tax_table slot or a list of two objects if behavior is "return_matrix"
#' @export
#' @examples
#' \dontrun{
#' # /!\ The value of threshold must be change for real database (recommend
#' #  value are between 50 and 70).
#'
#' data_fungi_mini_new <- assign_idtaxa(data_fungi_mini,
#'   fasta_for_training = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
#'     package = "MiscMetabar"
#'   ), threshold = 20, behavior = "add_to_phyloseq"
#' )
#'
#' result_idtaxa <- assign_idtaxa(data_fungi_mini,
#'   fasta_for_training = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
#'     package = "MiscMetabar"
#'   ), threshold = 20
#' )
#'
#' plot(result_idtaxa$idtaxa_raw)
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to [DECIPHER::IdTaxa()] if you
#'   use this function.
assign_idtaxa <- function(physeq,
                          trainingSet = NULL,
                          seq2search = NULL,
                          fasta_for_training = NULL,
                          behavior = "return_matrix",
                          threshold = 60,
                          column_names = c(
                            "Kingdom",
                            "Phylum",
                            "Class",
                            "Order",
                            "Family",
                            "Genus",
                            "Species"
                          ),
                          suffix = "_idtaxa",
                          nproc = 1,
                          unite = FALSE,
                          verbose = TRUE,
                          ...) {
  if (!is.null(trainingSet) && !is.null(fasta_for_training)) {
    stop("Please provide either trainingSet or fasta_for_training parameters, not both.")
  } else if (is.null(trainingSet) && is.null(fasta_for_training)) {
    stop("Please provide either trainingSet or fasta_for_training parameters.")
  } else if (is.null(trainingSet) && !is.null(fasta_for_training)) {
    if (verbose) {
      message("Training using fasta_for_training file.")
    }
    trainingSet_idtaxa <- learn_idtaxa(fasta_for_training, unite = unite)
  } else if (!is.null(trainingSet_idtaxa)) {
    trainingSet_idtaxa <- trainingSet
  }

  fasta2search <- write_temp_fasta(
    physeq = physeq,
    seq2search = seq2search,
    return_DNAStringSet = TRUE
  )

  fasta2search <- DECIPHER::OrientNucleotides(DECIPHER::RemoveGaps(fasta2search))

  if (verbose) {
    message("Classifing using training Set with IdTaxa.")
  }

  idtaxa_taxa_test <- DECIPHER::IdTaxa(
    test = fasta2search,
    trainingSet = trainingSet_idtaxa,
    processors = nproc,
    threshold = threshold,
    ...
  )

  idtaxa_taxa_df <- sapply(
    idtaxa_taxa_test,
    function(x) {
      paste(x$taxon,
        collapse = ";"
      )
    }
  )

  col2add <- max(stringr::str_count(idtaxa_taxa_df, ";")) - stringr::str_count(idtaxa_taxa_df, ";")
  for (i in seq_along(idtaxa_taxa_df)) {
    idtaxa_taxa_df[i] <-
      paste0(idtaxa_taxa_df[i], paste(as.character(rep(";", each = col2add[i])), collapse = ""))
  }
  t_idtaxa <- tibble::tibble(data.frame(stringr::str_split_fixed(idtaxa_taxa_df, ";", max(stringr::str_count(idtaxa_taxa_df, ";")) + 1)))[, -1]

  column_names <- paste0(column_names, suffix)
  colnames(t_idtaxa) <- column_names
  t_idtaxa$taxa_names <- names(fasta2search)

  if (behavior == "return_matrix") {
    return(list(
      "taxo_value" = t_idtaxa,
      "idtaxa_raw" = idtaxa_taxa_test
    ))
  } else if (behavior == "add_to_phyloseq") {
    tax_tab <- as.data.frame(as.matrix(physeq@tax_table))
    tax_tab$taxa_names <- taxa_names(physeq)

    new_physeq <- physeq

    new_tax_tab <- left_join(tax_tab, t_idtaxa,
      by = join_by(taxa_names)
    ) |>
      dplyr::select(-taxa_names) |>
      as.matrix()

    new_physeq@tax_table <- tax_table(new_tax_tab)
    taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

    return(new_physeq)
  } else {
    stop("Param behavior must take either 'return_matrix' or 'add_to_phyloseq' value")
  }
}
################################################################################

################################################################################
#' A wrapper of [DECIPHER::LearnTaxa()]
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' This function is basically a wrapper of functions [DECIPHER::LearnTaxa()],
#'  please cite the DECIPHER package if you use this function.
#'
#' @param fasta_for_training A fasta file (can be gzip) to train the trainingSet
#'   using the function [learn_idtaxa()]. Only used if trainingSet is NULL.
#'
#'   The reference database must contain
#'   taxonomic information in the header of each sequence in the form of a string
#'   starting with ";tax=" and followed by a comma-separated list of up to nine
#'   taxonomic identifiers.
#'
#'   The only exception is if `unite=TRUE`. In that case the UNITE taxonomy is
#'   automatically formatted.
#' @param output_Rdata A vector naming the path to an output Rdata file. If
#'   left to NULL, no Rdata file is written.
#' @param output_path_only (logical, default FALSE). If TRUE, the function
#'   return only the path to the output_Rdata file. Note that output_Rdata must
#'   be set.
#' @param unite (logical, default FALSE). If set to TRUE, the fasta_for_training
#'   file is formatted from UNITE format to sintax one, needed in
#'   fasta_for_training. Only used if trainingSet is NULL.
#' @param ... Additional arguments passed on to [DECIPHER::LearnTaxa()]
#'
#' @seealso [assign_idtaxa()]
#' @author Adrien Taudière
#' @return Either a Taxa Train object (see [DECIPHER::LearnTaxa()]) or, if
#'   output_path_only is TRUE, a vector indicating the path to the output
#'   training object.
#' @export
#' @examples
#' \dontrun{
#' training_mini_UNITE_fungi <-
#'   learn_idtaxa(fasta_for_training = system.file("extdata",
#'     "mini_UNITE_fungi.fasta.gz",
#'     package = "MiscMetabar"
#'   ))
#' plot(training_mini_UNITE_fungi)
#'
#' training_100sp_UNITE <-
#'   learn_idtaxa(
#'     fasta_for_training = system.file("extdata",
#'       "100_sp_UNITE_sh_general_release_dynamic.fasta",
#'       package = "MiscMetabar"
#'     ),
#'     unite = TRUE
#'   )
#'
#' plot(training_100sp_UNITE)
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to [DECIPHER::LearnTaxa()] if you
#'   use this function.
learn_idtaxa <- function(fasta_for_training, output_Rdata = NULL, output_path_only = FALSE, unite = FALSE, ...) {
  seqs <- Biostrings::readDNAStringSet(fasta_for_training)
  seqs <- DECIPHER::RemoveGaps(seqs)
  seqs <- DECIPHER::OrientNucleotides(seqs)

  taxo_for_learning <- names(seqs)
  if (unite) {
    taxo_for_learning <- gsub("(.*)(FU|reps)", "Root;", taxo_for_learning)
    taxo_for_learning <- gsub("(.*)(FU|reps_singleton)", "Root;", taxo_for_learning)
    taxo_for_learning <- gsub("(.*)(FU|refs_singleton)", "Root;", taxo_for_learning)
    taxo_for_learning <- gsub("(.*)(FU|refs)", "Root;", taxo_for_learning)
    taxo_for_learning <- gsub("|", "", taxo_for_learning, fixed = TRUE)
  }
  taxo_for_learning <-
    gsub(",", ";", gsub(
      ";tax=",
      "Root;",
      gsub("(.*)(;tax=)", "\\2", taxo_for_learning)
    )) # extract the group label

  taxo_for_learning[!grepl("^Root;", taxo_for_learning)] <-
    paste0("Root;", taxo_for_learning[!grepl("^Root;", taxo_for_learning)])

  train_idtaxa <- DECIPHER::LearnTaxa(seqs, taxonomy = taxo_for_learning, ...)
  if (!is.null(output_Rdata)) {
    save(train_idtaxa, output_Rdata)
  }
  if (output_path_only) {
    if (is.null(output_Rdata)) {
      stop("Param output_Rdata must be fill if output_path_only is TRUE.")
    }
    return(output_Rdata)
  } else {
    return(train_idtaxa)
  }
}
################################################################################


################################################################################

#' Assign taxonomy with dada2 using 2 steps assignTaxonomy and assignSpecies
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Mainly a wrapper of [dada2::assignTaxonomy()] and [dada2::assignSpecies()]
#'
#' @inheritParams clean_pq
#' @param ref_fasta (required) A link to a database in fasta.
#' @param seq2search A DNAStringSet object of sequences to search for. Replace
#'   the physeq object.
#' @param min_bootstrap (Float \[0:1\], default 0.5), See [dada2::assignTaxonomy()]
#' @param tryRC See [dada2::assignTaxonomy()]
#' @param taxa_ranks (vector of character) names for the column of the
#'   taxonomy
#' @param suffix (character) The suffix to name the new columns.
#'   Default to "_idtaxa".
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use
#' @param use_assignSpecies (logical, default TRUE) Do the Species rank is
#'   obtained using [dada2::assignSpecies()] ?
#' @param trunc_absent_ranks (logical, default FALSE) Do ranks present in
#'   taxa_ranks but not present in the database are removed ?
#' @param nproc (Float \[0:1\], default 0.5)
#' @param verbose (logical). If TRUE, print additional information.
#' @param seq_at_one_time How many sequences are treated at one time.
#'   See param `n` in [dada2::assignSpecies()]
#' @param allowMultiple (logical, default FALSE). Unchanged from
#'  [dada2::assignSpecies()].
#'  Defines the behavior when multiple exact matches
#'  against different species are returned. By default only unambiguous
#'  identifications are return. If TRUE, a concatenated string of all
#'  exactly matched species is returned. If an integer is provided,
#'  multiple identifications up to that many are returned as
#' a concatenated string.
#' @param from_sintax (logical, default FALSE). Set to TRUE
#'   if the ref_fasta database is in sintax format. See [assign_sintax()]
#'   for more information about the sintax format.
#'
#' @return Either a an object of class phyloseq (if `physeq` is not NULL),
#'   or a taxonomic table if `seq2search` is used in place of `physeq`
#' @export
#'
#' @examples
#' \dontrun{
#' data_fungi_mini2 <- assign_dada2(data_fungi_mini,
#'   ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz",
#'     package = "MiscMetabar"
#'   ), suffix = "_dada2",
#'   from_sintax = TRUE
#' )
#' }
assign_dada2 <- function(physeq = NULL,
                         ref_fasta = NULL,
                         seq2search = NULL,
                         min_bootstrap = 0.5,
                         tryRC = FALSE,
                         taxa_ranks = c(
                           "Kingdom",
                           "Phylum", "Class", "Order", "Family",
                           "Genus", "Species", "taxId"
                         ),
                         use_assignSpecies = TRUE,
                         trunc_absent_ranks = FALSE,
                         nproc = 1,
                         suffix = "",
                         verbose = TRUE,
                         seq_at_one_time = 2000,
                         allowMultiple = FALSE,
                         from_sintax = FALSE) {
  if (is.null(physeq)) {
    if (is.null(seq2search)) {
      stop("You must specify either physeq or seq2search args.")
    }
  } else {
    verify_pq(physeq)
    seq2search <- physeq@refseq
  }
  if (from_sintax) {
    Biostrings::writeXStringSet(format2dada2(fasta_db = ref_fasta), "temp_taxo.fasta")
    ref_fasta_taxo <- "temp_taxo.fasta"
  } else {
    ref_fasta_taxo <- ref_fasta
  }
  taxtab <- assignTaxonomy(
    seqs = seq2search,
    refFasta = ref_fasta_taxo,
    taxLevels = taxa_ranks,
    multithread = nproc,
    minBoot = min_bootstrap * 100,
    tryRC = tryRC,
    verbose = verbose,
    outputBootstraps = FALSE
  )

  taxtab <- as.data.frame(taxtab)
  rownames(taxtab) <- names(seq2search)

  if (use_assignSpecies) {
    Biostrings::writeXStringSet(format2dada2_species(fasta_db = ref_fasta), "temp_species.fasta")

    GS <- assignSpecies(
      seqs = seq2search,
      refFasta = "temp_species.fasta",
      allowMultiple = allowMultiple,
      n = seq_at_one_time,
      tryRC = tryRC
    )

    taxtab$Species <- as.vector(paste(GS[, 1], GS[, 2], sep = "_")) %>%
      {
        gsub("NA_NA", NA, .)
      }

    unlink("temp_species.fasta")
  }

  if (from_sintax) {
    unlink("temp_taxo.fasta")
  }

  if (!trunc_absent_ranks) {
    for (rank in taxa_ranks[!taxa_ranks %in% colnames(taxtab)]) {
      taxtab[[rank]] <- rep(NA, nrow(taxtab))
    }
  }

  colnames(taxtab) <- paste0(colnames(taxtab), suffix)
  if (is.null(physeq)) {
    return(taxtab)
  } else {
    tax_tab <- as.data.frame(as.matrix(physeq@tax_table))
    tax_tab$taxa_names <- taxa_names(physeq)

    taxtab$taxa_names <- rownames(taxtab)

    new_physeq <- physeq
    new_tax_tab <- left_join(tax_tab, taxtab,
      by = join_by(taxa_names)
    ) |>
      dplyr::select(-taxa_names) |>
      as.matrix()
    new_physeq@tax_table <- tax_table(new_tax_tab)
    taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

    return(new_physeq)
  }
}
################################################################################


################################################################################
#' Filter taxa by cleaning taxa with NA at given taxonomic rank(s)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Basically a wrapper of subset_taxa_pq()
#'
#' @inheritParams clean_pq
#' @param taxa_ranks A vector of taxonomic ranks. For examples c("Family","Genus").
#'   If taxa_ranks is NULL (default), all ranks are used, i.e. all taxa with at least 1 NA
#'   will be filtered out. Numeric position of taxonomic ranks can also be used.
#' @param n_NA (int default = 0). Number of allowed NA by taxa in the list of the
#'   taxonomic ranks
#' @param verbose (logical). If TRUE, print additional information.
#' @param NA_equivalent (vector of character, default NULL). Exact matching of
#'   the character listed in the vector are converted as NA before to filter out
#'   taxa.
#' @param clean_pq (logical, default TRUE)
#'   If set to TRUE, empty samples are discarded after filtering. See [clean_pq()].
#'
#' @return An object of class phyloseq
#'
#' @seealso [subset_taxa_pq()]
#' @export
#' @author Adrien Taudière
#'
#' @examples
#'
#' data_fungi_wo_NA <- filt_taxa_wo_NA(data_fungi)
#' filt_taxa_wo_NA(data_fungi, n_NA = 1)
#' filt_taxa_wo_NA(data_fungi, taxa_ranks = c(1:3))
#'
#' filt_taxa_wo_NA(data_fungi, taxa_ranks = c("Trait", "Confidence.Ranking"))
#' filt_taxa_wo_NA(data_fungi,
#'   taxa_ranks = c("Trait", "Confidence.Ranking"),
#'   NA_equivalent = c("-", "NULL")
#' )
filt_taxa_wo_NA <- function(physeq,
                            taxa_ranks = NULL,
                            n_NA = 0,
                            verbose = TRUE,
                            NA_equivalent = NULL,
                            clean_pq = TRUE) {
  verify_pq(physeq)
  new_physeq <- physeq

  if (is.null(taxa_ranks)) {
    taxa_ranks <- c(1:ncol(physeq@tax_table))
  }

  taxtab <- data.frame(unclass(physeq@tax_table[, taxa_ranks]))

  if (!is.null(NA_equivalent)) {
    taxtab <- taxtab |>
      mutate(across(everything(), ~ if_else(. %in% NA_equivalent, NA, .)))
  }

  cond <- rowSums(is.na(taxtab))
  new_physeq <- subset_taxa_pq(new_physeq, cond < (n_NA + 1), clean_pq = clean_pq)

  if (verbose) {
    message(
      "You filtered out ",
      ntaxa(physeq) - ntaxa(new_physeq),
      " taxa, leading to a phyloseq object including ",
      ntaxa(new_physeq),
      " taxa without NA in the taxonomic ranks: ",
      paste(taxa_ranks, collapse = " "),
      "."
    )
  }

  return(new_physeq)
}
################################################################################
