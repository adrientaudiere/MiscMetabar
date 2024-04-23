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
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - temp.fasta (refseq in fasta)
#'   - cluster.fasta (centroid)
#'   - temp.uc (clusters)
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
#' file_dna <- tempfile("dna.fa")
#' seqinr::write.fasta("GCCCATTAGTATTCTAGTGGGCATGCCTGTTCGAGCGTCATTTTCAACC",
#'   file = file_dna, names = "seq1"
#' )
#'
#' res <- vs_search_global(data_fungi, path_to_fasta = file_dna)
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
#'   Please cite [vsearch](https://github.com/torognes/vsearch).
#' @author Adrien Taudière

vs_search_global <- function(physeq,
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
    if (file.exists(paste0(tempdir(), "temp.fasta"))) {
      unlink(paste0(tempdir(), "temp.fasta"))
    }
    if (file.exists(paste0(tempdir(), "temp.uc"))) {
      unlink(paste0(tempdir(), "temp.uc"))
    }
    if (file.exists(paste0(tempdir(), "seq2search.fasta"))) {
      unlink(paste0(tempdir(), "seq2search.fasta"))
    }
  } else {
    message(paste0("Temporary files are located at ", tempdir()))
  }

  return(invisible(pack_clusts))
}
################################################################################



###############################################################################
#' Re-cluster sequences of an object of class `physeq`
#'   or cluster a list of DNA sequences using SWARM
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @inheritParams clean_pq
#' @param dna_seq NOT WORKING FOR THE MOMENT
#'   You may directly use a character vector of DNA sequences
#'   in place of physeq args. When physeq is set, dna sequences take the value of
#'   `physeq@refseq`
#' @param d (default: 1) maximum number of differences allowed between two
#'   amplicons, meaning that two amplicons will be grouped if they have `d`
#'   (or less) differences
#' @param swarmpath (default: swarm) path to swarm
#' @param vsearch_path (default: vsearch) path to vsearch, used only if physeq
#'   is NULL and dna_seq is provided.
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param swarm_args (default : "--fastidious") a one length character
#'   element defining other parameters to  passed on to swarm See other possible
#'   methods in the [SWARM pdf manual](https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf)
#' @param tax_adjust (Default 0) See the man page
#'   of [merge_taxa_vec()] for more details.
#'   To conserved the taxonomic rank of the most abundant ASV,
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary
#'   files ?
#'   - temp.fasta (refseq in fasta or dna_seq sequences)
#'   - temp_output (classical output of SWARM)
#'   - temp_uclust (clusters output of SWARM)
#' @details This function use the `merge_taxa_vec` function to
#'   merge taxa into clusters. By default tax_adjust = 0. See the man page
#'   of [merge_taxa_vec()].
#' @return A new object of class `physeq` or a list of cluster if dna_seq
#'   args was used.
#'
#' @references
#'   SWARM can be downloaded from
#'   \url{https://github.com/torognes/swarm/}.
#'
#' @export
#' @examplesIf MiscMetabar::is_swarm_installed()
#' summary_plot_pq(data_fungi)
#' system2("swarm", "-h")
#'
#' data_fungi_swarm <- swarm_clustering(data_fungi)
#' summary_plot_pq(data_fungi_swarm)
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
#'
#' sequences_ex_swarm <- swarm_clustering(
#'   dna_seq = sequences_ex
#' )
#' @seealso [asv2otu()], [vsearch_clustering()]
#' @references
#'   SWARM can be downloaded from
#'   \url{https://github.com/torognes/swarm}.
#'   More information in the associated publications
#'   \doi{doi:10.1093/bioinformatics/btab493} and \doi{doi:10.7717/peerj.593}
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [SWARM](https://github.com/torognes/swarm).

swarm_clustering <- function(physeq = NULL,
                             dna_seq = NULL,
                             d = 1,
                             swarmpath = "swarm",
                             vsearch_path = "vsearch",
                             nproc = 1,
                             swarm_args = "--fastidious",
                             tax_adjust = 0,
                             keep_temporary_files = FALSE) {
  dna <- physeq_or_string_to_dna(
    physeq = physeq,
    dna_seq = dna_seq
  )

  if (!is.null(physeq)) {
    nseq <- taxa_sums(physeq)
    nseq <- nseq[match(names(nseq), names(dna))]
    names(dna) <- paste0(names(dna), "_", nseq)
    Biostrings::writeXStringSet(dna, paste0(tempdir(), "/", "temp.fasta"))
    system2(
      swarmpath,
      paste0(
        paste0(tempdir(), "/", "temp.fasta"),
        " -o ",
        paste0(tempdir(), "/", "temp_output"),
        " ",
        " -u ",
        paste0(tempdir(), "/", "temp_uclust"),
        " -t ",
        nproc,
        " ",
        swarm_args
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  } else {
    Biostrings::writeXStringSet(dna, paste0(tempdir(), "/", "amplicons.fasta"))
    system2(
      vsearch_path,
      paste0(
        " --derep_fulllength ",
        paste0(tempdir(), "/", "amplicons.fasta"),
        " --sizeout  --relabel_sha1 --fasta_width 0 --output ",
        paste0(tempdir(), "/", "temp.fasta")
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    system2(
      swarmpath,
      paste0(
        paste0(tempdir(), "/", "temp.fasta"),
        " -o ",
        paste0(tempdir(), "/", "temp_output"),
        " ",
        " -u ",
        paste0(tempdir(), "/", "temp_uclust"),
        " -z ",
        " -t ",
        nproc,
        " ",
        swarm_args
      ),
      stdout = TRUE,
      stderr = TRUE
    )
  }

  pack_clusts <-
    utils::read.table(paste0(tempdir(), "/", "temp_uclust"), sep = "\t")
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

  if (inherits(physeq, "phyloseq")) {
    clusters <- pack_clusts$cluster[pack_clusts$type != "C"]
    names(clusters) <-
      sub("_.*$", "", pack_clusts$query[pack_clusts$type != "C"])

    clusters <- clusters[match(taxa_names(physeq), names(clusters))]

    new_obj <-
      merge_taxa_vec(physeq,
        clusters,
        tax_adjust = tax_adjust
      )
  } else if (inherits(dna_seq, "character")) {
    new_obj <- pack_clusts
  } else {
    stop(
      "You must set the args physeq (object of class phyloseq) or
      dna_seq (character vector)."
    )
  }

  if (file.exists(paste0(tempdir(), "/", "temp.fasta")) &&
    !keep_temporary_files) {
    unlink(paste0(tempdir(), "/", "temp.fasta"))
  }
  if (file.exists(paste0(tempdir(), "/", "temp_output")) &&
    !keep_temporary_files) {
    unlink(paste0(tempdir(), "/", "temp_output"))
  }
  if (file.exists(paste0(tempdir(), "/", "temp_uclust")) &&
    !keep_temporary_files) {
    unlink(paste0(tempdir(), "/", "temp_uclust"))
  }
  if (file.exists(paste0(tempdir(), "/", "amplicon.fasta")) &&
    !keep_temporary_files) {
    unlink(paste0(tempdir(), "/", "amplicon.fasta"))
  }
  return(new_obj)
}
###############################################################################



###############################################################################
#' Recluster sequences of an object of class `physeq`
#'   or cluster a list of DNA sequences using vsearch software
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @inheritParams clean_pq
#' @param dna_seq You may directly use a character vector of DNA sequences
#'   in place of physeq args. When physeq is set, dna sequences take the value of
#'   `physeq@refseq`
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param id (default: 0.97) level of identity to cluster
#' @param vsearchpath (default: vsearch) path to vsearch
#' @param tax_adjust (Default 0) See the man page
#'   of [merge_taxa_vec()] for more details.
#'   To conserved the taxonomic rank of the most abundant ASV,
#'   set tax_adjust to 0 (default). For the moment only tax_adjust = 0 is
#'   robust
#' @param vsearch_cluster_method (default: "--cluster_size) See other possible
#'   methods in the [vsearch manual](https://github.com/torognes/vsearch) (e.g. `--cluster_size` or `--cluster_smallmem`)
#'   - `--cluster_fast` : Clusterize the fasta sequences in filename, automatically sort by decreasing sequence length beforehand.
#'   - `--cluster_size` : Clusterize the fasta sequences in filename, automatically sort by decreasing sequence abundance beforehand.
#'   - `--cluster_smallmem` : Clusterize the fasta sequences in filename without automatically modifying their order beforehand. Sequence are expected to be sorted by decreasing sequence length, unless *--usersort* is used.
#'     In that case you may set `vsearch_args` to vsearch_args = "--strand both --usersort"
#' @param vsearch_args (default : "--strand both") a one length character element defining other parameters to
#'   passed on to vsearch.
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files ?
#'   - temp.fasta (refseq in fasta or dna_seq sequences)
#'   - cluster.fasta (centroid if method = "vsearch")
#'   - temp.uc (clusters if method = "vsearch")
#'
#' @seealso [asv2otu()], [swarm_clustering()]
#' @details This function use the [merge_taxa_vec()] function to
#'   merge taxa into clusters. By default tax_adjust = 0. See the man page
#'   of [merge_taxa_vec()].
#'
#' @return A new object of class `physeq` or a list of cluster if dna_seq
#'   args was used.
#'
#' @references
#'   VSEARCH can be downloaded from
#'   \url{https://github.com/torognes/vsearch}.
#'   More information in the associated publication
#'   \url{https://pubmed.ncbi.nlm.nih.gov/27781170}.
#' @export
#' @author Adrien Taudière
#'
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
#' summary_plot_pq(data_fungi)
#' d_vs <- vsearch_clustering(data_fungi)
#' summary_plot_pq(d_vs)
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [vsearch](https://github.com/torognes/vsearch).
vsearch_clustering <- function(physeq = NULL,
                               dna_seq = NULL,
                               nproc = 1,
                               id = 0.97,
                               vsearchpath = "vsearch",
                               tax_adjust = 0,
                               vsearch_cluster_method = "--cluster_size",
                               vsearch_args = "--strand both",
                               keep_temporary_files = FALSE) {
  dna <- physeq_or_string_to_dna(physeq = physeq, dna_seq = dna_seq)

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
  } else if (inherits(dna_seq, "character")) {
    new_obj <- pack_clusts
  } else {
    stop(
      "You must set the args physeq (object of class phyloseq) or dna_seq (character vector)."
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
  return(new_obj)
}
###############################################################################


################################################################################
#' Search for a list of sequence in an object to remove chimera taxa
#'   using [vsearch](https://github.com/torognes/vsearch)
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param object (required) A phyloseq-class object or one of dada, derep,
#'   data.frame or list coercible to sequences table using the
#'   function [dada2::makeSequenceTable()]
#' @param type (default "Discard_only_chim"). The type define the type of
#'   filtering.
#'
#'  - "Discard_only_chim" will only discard taxa classify as chimera by vsearch
#'  - "Select_only_non_chim" will only select taxa classify as non-chimera by
#'    vsearch(after filtering taxa based on their sequence length by the
#'    parameter `min_seq_length` from the [chimera_detection_vs()] function)
#'  - "Select_only_chim" will only select taxa classify as chimera by
#'    vsearch (after filtering taxa based on their sequence length by the
#'    parameter `min_seq_length` from the [chimera_detection_vs()] function)
#' @param clean_pq (logical; default FALSE) If TRUE, return the phyloseq object
#'   after cleaning using the default parameter of [clean_pq()] function.
#'
#' @param ... Others arguments passed on to [chimera_detection_vs()] function
#' @seealso [chimera_detection_vs()]
#' @return
#'
#' - I/ a sequences tables if object is of class dada, derep, data.frame or
#'   list.
#' - II/ a phyloseq object without (or with if type = 'Select_only_chim')
#'   chimeric taxa
#'
#' @export
#'
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
#' data_fungi_nochim <- chimera_removal_vs(data_fungi)
#' data_fungi_nochim_16 <- chimera_removal_vs(data_fungi,
#'   abskew = 16,
#'   min_seq_length = 10
#' )
#' data_fungi_nochim2 <-
#'   chimera_removal_vs(data_fungi, type = "Select_only_non_chim")
#' data_fungi_chimera <-
#'   chimera_removal_vs(data_fungi, type = "Select_only_chim")
#' }
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make [vsearch](https://github.com/torognes/vsearch).


chimera_removal_vs <-
  function(object,
           type = "Discard_only_chim",
           clean_pq = FALSE,
           ...) {
    if (inherits(object, "dada") ||
      inherits(object, "derep") ||
      inherits(object, "data.frame") ||
      inherits(object, "list")) {
      object <- makeSequenceTable(object)
    }

    if (inherits(object, "matrix")) {
      chim_detect <-
        chimera_detection_vs(
          seq2search = colnames(object),
          nb_seq = colSums(object),
          ...
        )
      if (type == "Discard_only_chim") {
        seq_tab_final <-
          object[, !colnames(object) %in% as.character(chim_detect$chimera)]
      } else if (type == "Select_only_non_chim") {
        seq_tab_final <- seq_tab_Pairs[, as.character(chim_rm$non_chimera)]
      } else if (type == "Select_only_chim") {
        seq_tab_final <- seq_tab_Pairs[, as.character(chim_rm$chimera)]
      } else {
        stop(
          "Type must be set to one of 'Discard_only_chim',
             'Select_only_non_chim', or 'Select_only_chim'"
        )
      }
      seq_tab_final <- seq_tab_final
      return(seq_tab_final)
    } else if (inherits(object, "phyloseq")) {
      verify_pq(object)
      if (sum(taxa_sums(object) == 0) > 0) {
        object <- clean_pq(object)
      }
      chim_detect <-
        chimera_detection_vs(
          seq2search = refseq(object),
          nb_seq = taxa_sums(object),
          ...
        )
      if (type == "Discard_only_chim") {
        cond <-
          !as.character(refseq(object)) %in% as.character(chim_detect$chimera)
        names(cond) <- names(refseq(object))
        new_physeq <-
          subset_taxa_pq(object, condition = cond)
      } else if (type == "Select_only_non_chim") {
        cond <-
          as.character(refseq(object)) %in% as.character(chim_detect$non_chimera)
        names(cond) <- names(refseq(object))
        new_physeq <-
          subset_taxa_pq(object, condition = cond)
      } else if (type == "Select_only_chim") {
        cond <-
          as.character(refseq(object)) %in% as.character(chim_detect$chimera)
        names(cond) <- names(refseq(object))
        new_physeq <-
          subset_taxa_pq(object, condition = cond)
      } else {
        stop(
          "Type must be set to one of 'Discard_only_chim',
             'Select_only_non_chim', or 'Select_only_chim'"
        )
      }
      if (clean_pq) {
        new_physeq <- clean_pq(physeq)
      }
      return(new_physeq)
    }
  }

#' Detect for chimera taxa using [vsearch](https://github.com/torognes/vsearch)
#'
#' `r lifecycle::badge("experimental")`
#'
#'
#' @param seq2search (required) a list of DNA sequences coercible by function
#'  [Biostrings::DNAStringSet()]
#' @param nb_seq (required) a numeric vector giving the number of sequences for
#'   each DNA sequences
#' @param vsearchpath (default: vsearch) path to vsearch
#' @param abskew (int, default 2) The abundance skew is used to distinguish in a
#'   three way alignment which sequence is the chimera and which are the
#'   parents. The assumption is that chimeras appear later in the PCR
#'   amplification process and are therefore less abundant than their parents.
#'   The default value is 2.0, which means that the parents should be at least
#'   2 times more abundant than their chimera. Any positive value equal or
#'   greater than 1.0 can be used.
#' @param min_seq_length (int, default 100)) Minimum length of sequences to
#'   be part of the analysis
#' @param vsearch_args (default "--fasta_width 0") A list of other args for
#'   vsearch command
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary
#'  files ?
#'   - non_chimeras.fasta
#'   - chimeras.fasta
#'   - borderline.fasta
#'
#' @return A list of 3 including non-chimera taxa (`$non_chimera`), chimera taxa
#'   (`$chimera`) and bordeline taxa (`$borderline`)
#' @export
#'
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
#' chimera_detection_vs(
#'   seq2search = data_fungi@refseq,
#'   nb_seq = taxa_sums(data_fungi)
#' )
#' }
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make [vsearch](https://github.com/torognes/vsearch).
chimera_detection_vs <- function(seq2search,
                                 nb_seq,
                                 vsearchpath = "vsearch",
                                 abskew = 2,
                                 min_seq_length = 100,
                                 vsearch_args = "--fasta_width 0",
                                 keep_temporary_files = FALSE) {
  dna_raw <- Biostrings::DNAStringSet(seq2search)
  names(dna_raw) <- paste0(
    "ASV", seq(1, length(seq2search)),
    ";size=", nb_seq
  )

  dna <- dna_raw[Biostrings::width(dna_raw) >= min_seq_length]
  abun <- unlist(strsplit(names(dna), split = "="))
  abun_tot <-
    sum(as.numeric(abun[seq(2, 2 * length(abun), by = 2)]), na.rm = T)

  message(
    paste(
      "Filtering for sequences under",
      min_seq_length,
      "bp remove a total of",
      length(dna_raw) - length(dna),
      "(",
      round((length(dna_raw) - length(dna)) / length(dna_raw) * 100, 2),
      "%)",
      "unique sequences for a total of",
      sum(nb_seq) - abun_tot,
      "sequences removed",
      "(",
      round((sum(nb_seq) - abun_tot) / sum(nb_seq) * 100, 2),
      "%)"
    )
  )

  Biostrings::writeXStringSet(
    dna,
    paste0(tempdir(), "/", "temp.fasta")
  )

  system2(
    vsearchpath,
    paste0(
      " --uchime_denovo ",
      paste0(tempdir(), "/", "temp.fasta"),
      " --abskew ",
      abskew,
      " --nonchimeras ",
      paste0(tempdir(), "/", "non_chimeras.fasta"),
      " --chimeras ",
      paste0(tempdir(), "/", "chimeras.fasta"),
      " --borderline ",
      paste0(tempdir(), "/", "borderline.fasta"),
      " ",
      vsearch_args
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  non_chimera_AAStringSet <-
    Biostrings::readAAStringSet(paste0(tempdir(), "/", "non_chimeras.fasta"))
  chimera_AAStringSet <-
    Biostrings::readAAStringSet(paste0(tempdir(), "/", "chimeras.fasta"))
  borderline_AAStringSet <-
    Biostrings::readAAStringSet(paste0(tempdir(), "/", "borderline.fasta"))

  if (!keep_temporary_files) {
    if (file.exists(paste0(tempdir(), "temp.fasta"))) {
      unlink(paste0(tempdir(), "temp.fasta"))
    }
    if (file.exists(paste0(tempdir(), "non_chimeras.fasta"))) {
      unlink(paste0(tempdir(), "non_chimeras.fasta"))
    }
    if (file.exists(paste0(tempdir(), "chimeras.fasta"))) {
      unlink(paste0(tempdir(), "chimeras.fasta"))
    }
  } else {
    message(paste0("Temporary files are located at ", tempdir()))
  }

  return(
    list(
      "non_chimera" = non_chimera_AAStringSet,
      "chimera" = chimera_AAStringSet,
      "borderline" = borderline_AAStringSet
    )
  )
}
