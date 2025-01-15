################################################################################
#' Search for a list of sequence in a fasta file against physeq reference
#'   sequences using [vsearch](https://github.com/torognes/vsearch)
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Use of VSEARCH software.
#'
#' @inheritParams clean_pq
#' @param seq2search (required if path_to_fasta is NULL) Either (i) a DNAstringSet object
#'   or (ii) a character vector that will be convert to DNAstringSet using
#'   [Biostrings::DNAStringSet()]
#' @param path_to_fasta (required if seq2search is NULL) a path to fasta file if seq2search is est to NULL.
#' @param vsearchpath (default: "vsearch") path to vsearch
#' @param id (default: 0.8) id for the option `--usearch_global` of the vsearch software
#' @param iddef (default: 0) iddef for the option `--usearch_global` of the vsearch software
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files
#'   - temp.fasta (refseq in fasta)
#'   - cluster.fasta (centroid)
#'   - temp.uc (clusters)
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
#' if (requireNamespace("seqinr")) {
#'   file_dna <- tempfile("dna.fa")
#'   seqinr::write.fasta("GCCCATTAGTATTCTAGTGGGCATGCCTGTTCGAGCGTCATTTTCAACC",
#'     file = file_dna, names = "seq1"
#'   )
#'
#'   res <- vs_search_global(data_fungi, path_to_fasta = file_dna)
#'   unlink(file_dna)
#'
#'   res[res$identity != "*", ]
#'
#'   clean_pq(subset_taxa(data_fungi, res$identity != "*"))
#' }
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

  if (keep_temporary_files) {
    message(paste0("Temporary files are located at ", tempdir()))
  } else {
    if (file.exists(paste0(tempdir(), "temp.fasta"))) {
      unlink(paste0(tempdir(), "temp.fasta"))
    }
    if (file.exists(paste0(tempdir(), "temp.uc"))) {
      unlink(paste0(tempdir(), "temp.uc"))
    }
    if (file.exists(paste0(tempdir(), "seq2search.fasta"))) {
      unlink(paste0(tempdir(), "seq2search.fasta"))
    }
  }

  return(invisible(pack_clusts))
}
################################################################################



###############################################################################
#' Re-cluster sequences of an object of class `physeq`
#'   or cluster a list of DNA sequences using SWARM
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' A wrapper of SWARM software.
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
#' @seealso [postcluster_pq()], [vsearch_clustering()]
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' A wrapper of VSEARCH software.
#'
#' @inheritParams clean_pq
#' @param dna_seq You may directly use a character vector of DNA sequences
#'   in place of physeq args. When physeq is set, dna sequences take the value of
#'   `physeq@refseq`
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param id (default: 0.97) level of identity to cluster
#' @param vsearchpath (default: "vsearch") path to vsearch
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
#' @seealso [postcluster_pq()], [swarm_clustering()]
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
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use the VSEARCH software.
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
#' @param ... Other arguments passed on to [chimera_detection_vs()] function
#' @seealso [chimera_detection_vs()], [dada2::removeBimeraDenovo()]
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
################################################################################



################################################################################
#' Detect for chimera taxa using [vsearch](https://github.com/torognes/vsearch)
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Use the VSEARCH software.
#'
#' @param seq2search (required) a list of DNA sequences coercible by function
#'  [Biostrings::DNAStringSet()]
#' @param nb_seq (required) a numeric vector giving the number of sequences for
#'   each DNA sequences
#' @param vsearchpath (default: "vsearch") path to vsearch
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
#' @seealso [chimera_removal_vs()], [dada2::removeBimeraDenovo()]
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
    "Taxa", seq(1, length(seq2search)),
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

  if (keep_temporary_files) {
    message(paste0("Temporary files are located at ", tempdir()))
  } else {
    if (file.exists(paste0(tempdir(), "temp.fasta"))) {
      unlink(paste0(tempdir(), "temp.fasta"))
    }
    if (file.exists(paste0(tempdir(), "non_chimeras.fasta"))) {
      unlink(paste0(tempdir(), "non_chimeras.fasta"))
    }
    if (file.exists(paste0(tempdir(), "chimeras.fasta"))) {
      unlink(paste0(tempdir(), "chimeras.fasta"))
    }
  }

  return(
    list(
      "non_chimera" = non_chimera_AAStringSet,
      "chimera" = chimera_AAStringSet,
      "borderline" = borderline_AAStringSet
    )
  )
}
################################################################################

################################################################################
#' Write a temporary fasta file (internal use)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Write a fasta file from either a biostring object seq2search or a
#' refseq slot from a phyloseq object
#'
#' @inheritParams assign_sintax
#' @param temporary_fasta_file The name of a temporary_fasta_file (default "temp.fasta")
#' @param return_DNAStringSet (Logical default FALSE). If true, the temporary fasta file
#'   is removed and a DNAStringSet is return
#' @seealso [assign_sintax()], [assign_vsearch_lca]
#' @return Nothing, produce a fasta file or return a DNAStringset if temporary_fasta_file
#' @keywords internal
#' @noRd
#'

write_temp_fasta <- function(physeq,
                             seq2search,
                             temporary_fasta_file = "temp.fasta",
                             behavior = NULL,
                             clean_pq = TRUE,
                             verbose = TRUE,
                             return_DNAStringSet = FALSE) {
  if (!is.null(physeq) && !is.null(seq2search)) {
    stop("You must enter a single parameter from physeq and seq2search.")
  } else if (is.null(seq2search)) {
    verify_pq(physeq)
    if (is.null(physeq@refseq)) {
      stop("The phyloseq object do not contain a @refseq slot")
    }
    if (clean_pq) {
      physeq <- clean_pq(physeq, silent = !verbose)
    }
    dna <- Biostrings::DNAStringSet(physeq@refseq)
    Biostrings::writeXStringSet(dna, temporary_fasta_file)
  } else if (is.null(physeq)) {
    Biostrings::writeXStringSet(seq2search, temporary_fasta_file)
    if (behavior == "add_to_phyloseq") {
      stop("You can't use behavior = 'add_to_phyloseq' with seq2search param.")
    }
  } else if (is.null(physeq) && is.null(seq2search)) {
    stop("You must specify either physeq or seq2search parameter.")
  }
  if (return_DNAStringSet) {
    res <- Biostrings::readDNAStringSet(temporary_fasta_file)
    unlink(temporary_fasta_file)
    return(res)
  }
}
################################################################################

################################################################################
#' Assign Taxonomy using Sintax algorithm of Vsearch
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Please cite [Vsearch](https://github.com/torognes/vsearch)
#'   if you use this function to assign taxonomy.
#'
#' @inheritParams clean_pq
#' @param seq2search A DNAStringSet object of sequences to search for.
#' @param ref_fasta (required) A link to a database in vsearch format
#'  The reference database must contain taxonomic information in the header of
#'  each sequence in the form of a string starting with ";tax=" and followed
#'  by a comma-separated list of up to nine taxonomic identifiers. Each taxonomic
#'  identifier must start with an indication of the rank by one of the letters d
#'  (for domain) k (kingdom), p (phylum), c (class), o (order), f (family),
#'   g (genus), s (species), or t (strain). The letter is followed by a colon
#'    (:) and the name of that rank. Commas and semicolons are not allowed in
#'    the name of the rank. Non-ascii characters should be avoided in the names.
#'
#'  Example:
#'
#'  \>X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli,t:str._K-12_substr._MG1655
#' @param behavior Either "return_matrix" (default), "return_cmd",
#' or "add_to_phyloseq":
#'
#'  - "return_matrix" return a list of two matrix with taxonomic value in the
#'    first element of the list and bootstrap value in the second one.
#'
#'  - "return_cmd" return the command to run without running it.
#'
#'  - "add_to_phyloseq" return a phyloseq object with amended slot `@taxtable`.
#'    Only available if using physeq input and not seq2search input.
#'
#' @param vsearchpath (default: "vsearch") path to vsearch
#' @param clean_pq (logical, default TRUE)
#'   If set to TRUE, empty samples and empty ASV are discarded
#'   before clustering.
#' @param nproc (default: 1)
#'   Set to number of cpus/processors to use
#' @param suffix (character) The suffix to name the new columns.
#'   If set to "" (the default), the taxo_rank algorithm is used
#'   without suffix.
#' @param taxo_rank A list with the name of the taxonomic rank present in
#'   ref_fasta
#' @param min_boostrap (Int. \[0:1\], default 0.5)
#'   Minimum bootstrap value to inform taxonomy. For each bootstrap
#'   below the min_boostrap value, the taxonomy information is set to NA.
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files?
#'
#' - temporary_fasta_file (default "temp.fasta") : the fasta file from physeq
#'   or seq2search
#'
#' - "output_taxo_vs.txt" : see Vsearch Manual for parameter --tabbedout
#'
#' @param verbose (logical). If TRUE, print additional information.
#' @param temporary_fasta_file The name of a temporary_fasta_file (default "temp.fasta")
#' @param cmd_args Other arguments to be passed on to vsearch sintax cmd.
#'   By default cmd_args is equal to "--sintax_random" as recommended by
#'   [Torognes](https://github.com/torognes/vsearch/issues/535).
#' @param too_few (default value "align_start") see [tidyr::separate_wider_delim()]
#' @param too_many (default value "drop") see [tidyr::separate_wider_delim()]
#' @return See param behavior
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
#' assign_sintax(data_fungi_mini,
#'   ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
#'   behavior = "return_cmd"
#' )
#'
#' data_fungi_mini_new <- assign_sintax(data_fungi_mini,
#'   ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
#'   behavior = "add_to_phyloseq"
#' )
#'
#' assignation_results <- assign_sintax(data_fungi_mini,
#'   ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar")
#' )
#'
#' left_join(
#'   tidyr::pivot_longer(assignation_results$taxo_value, -taxa_names),
#'   tidyr::pivot_longer(assignation_results$taxo_boostrap, -taxa_names),
#'   by = join_by(taxa_names, name),
#'   suffix = c("rank", "bootstrap")
#' ) |>
#'   mutate(name = factor(name,
#'     levels = c(
#'       "Kingdom", "Phylum", "Class",
#'       "Order", "Family", "Genus", "Species"
#'     )
#'   )) |>
#'   # mutate(valuerank = forcats::fct_reorder(valuerank,
#'   #   as.integer(name), .desc = TRUE)) |>
#'   ggplot(aes(valuebootstrap,
#'     valuerank,
#'     fill = name
#'   )) +
#'   geom_jitter(alpha = 0.8, aes(color = name)) +
#'   geom_boxplot(alpha = 0.3)
#' }
#' @export
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [vsearch](https://github.com/torognes/vsearch).
assign_sintax <- function(physeq = NULL,
                          seq2search = NULL,
                          ref_fasta = NULL,
                          behavior = c("return_matrix", "add_to_phyloseq", "return_cmd"),
                          vsearchpath = "vsearch",
                          clean_pq = TRUE,
                          nproc = 1,
                          suffix = "",
                          taxo_rank = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                          min_boostrap = 0.5,
                          keep_temporary_files = FALSE,
                          verbose = TRUE,
                          temporary_fasta_file = "temp.fasta",
                          cmd_args = "--sintax_random",
                          too_few = "align_start",
                          too_many = "drop") {
  behavior <- match.arg(behavior)

  write_temp_fasta(
    physeq = physeq,
    seq2search = seq2search,
    temporary_fasta_file = temporary_fasta_file,
    behavior = behavior,
    clean_pq = clean_pq,
    verbose = verbose
  )


  if (verbose) {
    message("Start Vsearch sintax")
  }

  cmd_sintax <-
    paste0(
      " --sintax ",
      temporary_fasta_file,
      " --db ",
      ref_fasta,
      " --tabbedout output_taxo_vs.txt ",
      " --threads ",
      nproc,
      " ",
      cmd_args
    )

  if (behavior == "return_cmd") {
    if (!keep_temporary_files) {
      unlink(temporary_fasta_file)
    }
    return("sintax" = paste0(vsearchpath, " ", cmd_sintax))
  }

  system2(vsearchpath,
    args = cmd_sintax,
    stdout = TRUE,
    stderr = TRUE
  )

  if (!file.exists("output_taxo_vs.txt")) {
    warning("No taxonomic assignation were maded.")
    if (!keep_temporary_files) {
      unlink(temporary_fasta_file)
    }
    if (behavior == "add_to_phyloseq") {
      return(physeq)
    } else {
      return(NULL)
    }
  }
  res_sintax <- read.csv("output_taxo_vs.txt", sep = "\t", header = F)
  taxa_names <- res_sintax$V1
  res_sintax <- tibble(res_sintax$V2, taxa_names)
  res_sintax <- res_sintax |>
    tidyr::separate_wider_delim(-taxa_names, names = paste0(taxo_rank, suffix), delim = ",", too_few = too_few) |>
    tidyr::pivot_longer(-taxa_names) |>
    tidyr::separate_wider_delim(
      value,
      names_sep = "",
      names = c("", "_bootstrap"),
      delim = "(",
      too_many = too_many
    ) |>
    mutate(across(value_bootstrap, ~ as.numeric(gsub(")", "", .x)))) |>
    mutate(across(value, ~ gsub(".:", "", .x)))

  res_sintax_wide_bootstrap <-
    res_sintax |>
    select(-value) |>
    tidyr::pivot_wider(names_from = name, values_from = value_bootstrap)

  res_sintax_wide_taxo <-
    res_sintax |>
    select(-value_bootstrap) |>
    tidyr::pivot_wider(names_from = name, values_from = value)

  if (!is.null(min_boostrap)) {
    res_sintax_wide_taxo_filter <- res_sintax_wide_taxo
    res_sintax_wide_taxo_filter[res_sintax_wide_bootstrap < min_boostrap] <- NA
  }

  if (!keep_temporary_files) {
    unlink(temporary_fasta_file)
    unlink("output_taxo_vs.txt")
  }

  if (behavior == "add_to_phyloseq") {
    tax_tab <- as.data.frame(as.matrix(physeq@tax_table))
    tax_tab$taxa_names <- taxa_names(physeq)

    new_physeq <- physeq
    new_tax_tab <- left_join(tax_tab, res_sintax_wide_taxo_filter,
      by = join_by(taxa_names)
    ) |>
      dplyr::select(-taxa_names) |>
      as.matrix()
    new_physeq@tax_table <- tax_table(new_tax_tab)
    taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

    return(new_physeq)
  } else if (behavior == "return_matrix") {
    return(list(
      "taxo_value" = res_sintax_wide_taxo,
      "taxo_boostrap" = res_sintax_wide_bootstrap
    ))
  }
}
################################################################################

################################################################################
#' Assign taxonomy using LCA *à la* [stampa](https://github.com/frederic-mahe/stampa)
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'  Please cite [Vsearch](https://github.com/torognes/vsearch) and
#'   [stampa](https://github.com/frederic-mahe/stampa) if you use this function
#'   to assign taxonomy.
#'
#'
#' @inheritParams clean_pq
#' @param seq2search A DNAStringSet object of sequences to search for.
#' @param ref_fasta (required) A link to a database in vsearch format
#'  The reference database must contain taxonomic information in the header of
#'  each sequence in the form of a string starting with ";tax=" and followed
#'  by a comma-separated list of up to nine taxonomic identifiers. Each taxonomic
#'  identifier must start with an indication of the rank by one of the letters d
#'  (for domain) k (kingdom), p (phylum), c (class), o (order), f (family),
#'   g (genus), s (species), or t (strain). The letter is followed by a colon
#'    (:) and the name of that rank. Commas and semicolons are not allowed in
#'    the name of the rank. Non-ascii characters should be avoided in the names.
#'
#'  Example:
#'
#'  \>X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli,t:str._K-12_substr._MG1655
#' @param behavior Either "return_matrix" (default), "return_cmd",
#' or "add_to_phyloseq":
#'
#'  - "return_matrix" return a list of two matrix with taxonomic value in the
#'    first element of the list and bootstrap value in the second one.
#'
#'  - "return_cmd" return the command to run without running it.
#'
#'  - "add_to_phyloseq" return a phyloseq object with amended slot `@taxtable`.
#'    Only available if using physeq input and not seq2search input.
#'
#' @param vsearchpath (default: "vsearch") path to vsearch
#' @param clean_pq (logical, default TRUE)
#'   If set to TRUE, empty samples and empty ASV are discarded
#'   before clustering.
#' @param taxo_rank A list with the name of the taxonomic rank present in
#'   ref_fasta
#' @param nproc (int, default: 1)
#'   Set to number of cpus/processors to use
#' @param suffix (character) The suffix to name the new columns.
#'   If set to "" (the default), the taxo_rank algorithm is used
#'   without suffix.
#' @param id (Int. \[0:1\] default 0.5). Default value is based on
#'   [stampa](https://github.com/frederic-mahe/stampa).
#'   See Vsearch Manual for parameter `--id`
#' @param lca_cutoff (int, default 1). Fraction of matching hits
#'  required for the last common ancestor (LCA) output. For example, a value
#'  of 0.9 imply that if less than 10% of assigned species are not congruent
#'  the taxonomy is filled.
#'  Default value is based on [stampa](https://github.com/frederic-mahe/stampa).
#'  See Vsearch Manual for parameter `--lca_cutoff`
#'
#'  Text from vsearch manual :
#'  "Adjust the fraction of matching hits required for the last
#'  common ancestor (LCA) output with the --lcaout option during searches.
#'  The default value is 1.0 which requires all hits to match at each taxonomic
#'  rank for that rank to be included. If a lower cutoff value is used,
#'  e.g. 0.95, a small fraction of non-matching hits are allowed while that
#'  rank will still be reported. The argument to this option must be larger
#'  than 0.5, but not larger than 1.0"
#' @param maxrejects (int, default: 32)
#'   Maximum number of non-matching target sequences to consider before
#'   stopping the search for a given query.
#'   Default value is based on [stampa](https://github.com/frederic-mahe/stampa)
#'   See Vsearch Manual for parameter `--maxrejects`.

#' @param top_hits_only (Logical, default TRUE)
#'  Only the top hits with an equally high percentage of identity between the query and
#'  database sequence sets are written to the output. If you set top_hits_only
#'  you may need to set a lower `maxaccepts` and/or `lca_cutoof`.
#'   Default value is based on [stampa](https://github.com/frederic-mahe/stampa)
#'   See Vsearch Manual for parameter `--top_hits_only`
#' @param maxaccepts (int, default: 0)
#'   Default value is based on [stampa](https://github.com/frederic-mahe/stampa).
#'   Maximum number of matching target sequences to accept before stopping the search
#'   for a given query.
#'   See Vsearch Manual for parameter `--maxaccepts`
#' @param keep_temporary_files (logical, default: FALSE) Do we keep temporary files?
#'
#' - temporary_fasta_file (default "temp.fasta") : the fasta file from physeq or
#'   seq2search
#'
#' - "out_lca.txt" : see Vsearch Manual for parameter --lcaout
#'
#' - "userout.txt" : see Vsearch Manual for parameter --userout
#'
#' @param verbose (logical). If TRUE, print additional information.
#' @param temporary_fasta_file Name of the temporary fasta file. Only useful
#'   with keep_temporary_files = TRUE.
#' @param cmd_args Other arguments to be passed on to vsearch usearch_global cmd.
#' @param too_few (default value "align_start") see [tidyr::separate_wider_delim()]
#' @return See param behavior
#' @seealso [assign_sintax()], [add_new_taxonomy_pq()]
#' @examplesIf MiscMetabar::is_vsearch_installed()
#' \donttest{
#' data_fungi_mini_new <- assign_vsearch_lca(data_fungi_mini,
#'   ref_fasta = system.file("extdata", "mini_UNITE_fungi.fasta.gz", package = "MiscMetabar"),
#'   lca_cutoff = 0.9
#' )
#' }
#' @export
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite [vsearch](https://github.com/torognes/vsearch) and
#'   [stampa](https://github.com/frederic-mahe/stampa)
assign_vsearch_lca <- function(physeq = NULL,
                               seq2search = NULL,
                               ref_fasta = NULL,
                               behavior = c("return_matrix", "add_to_phyloseq", "return_cmd"),
                               vsearchpath = "vsearch",
                               clean_pq = TRUE,
                               taxo_rank = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                               nproc = 1,
                               suffix = "",
                               id = 0.5,
                               lca_cutoff = 1,
                               maxrejects = 32,
                               top_hits_only = TRUE,
                               maxaccepts = 0,
                               keep_temporary_files = FALSE,
                               verbose = TRUE,
                               temporary_fasta_file = "temp.fasta",
                               cmd_args = "",
                               too_few = "align_start") {
  behavior <- match.arg(behavior)
  write_temp_fasta(
    physeq = physeq,
    seq2search = seq2search,
    temporary_fasta_file = temporary_fasta_file,
    behavior = behavior,
    clean_pq = clean_pq,
    verbose = verbose
  )

  cmd_usearch <-
    paste0(
      " --usearch_global temp.fasta --db ",
      ref_fasta,
      " --lcaout out_lca.txt -id ",
      id,
      " --threads ",
      nproc,
      " --userfields query+id+target",
      " --maxaccepts ",
      maxaccepts,
      " --maxrejects ",
      maxrejects,
      " --lca_cutoff  ",
      lca_cutoff,
      " --userout userout.txt ",
      cmd_args
    )

  if (top_hits_only) {
    cmd_usearch <-
      paste0(cmd_usearch, " --top_hits_only")
  }

  if (behavior == "return_cmd") {
    if (!keep_temporary_files) {
      unlink(temporary_fasta_file)
    }
    return("sintax" = paste0(vsearchpath, " ", cmd_usearch))
  }

  system2(vsearchpath,
    args = cmd_usearch,
    stdout = TRUE,
    stderr = TRUE
  )

  res_usearch <- read.csv("out_lca.txt", sep = "\t", header = F)
  taxa_names <- res_usearch$V1
  res_usearch <- tibble(res_usearch$V2, taxa_names)
  if (sum(is.na(res_usearch)) == nrow(res_usearch)) {
    message("None match were found using usearch global and id=", id)
    if (behavior == "add_to_phyloseq") {
      return(physeq)
    } else {
      return(NULL)
    }
  }

  res_usearch <- res_usearch |>
    tidyr::separate_wider_delim(-taxa_names,
      names = paste0(taxo_rank, suffix),
      delim = ",",
      too_few = too_few
    ) |>
    tidyr::pivot_longer(-taxa_names) |>
    mutate(across(value, ~ gsub(".:", "", .x)))

  res_usearch_wide_taxo <-
    res_usearch |>
    tidyr::pivot_wider(names_from = name, values_from = value)

  if (!keep_temporary_files) {
    unlink(temporary_fasta_file)
    unlink("out_lca.txt")
    unlink("userout.txt")
  }

  if (behavior == "add_to_phyloseq") {
    tax_tab <- as.data.frame(as.matrix(physeq@tax_table))
    tax_tab$taxa_names <- taxa_names(physeq)

    new_physeq <- physeq
    new_tax_tab <- left_join(tax_tab, res_usearch_wide_taxo,
      by = join_by(taxa_names)
    ) |>
      dplyr::select(-taxa_names) |>
      as.matrix()
    new_physeq@tax_table <- tax_table(new_tax_tab)
    taxa_names(new_physeq@tax_table) <- taxa_names(physeq)

    return(new_physeq)
  } else if (behavior == "return_matrix") {
    return(res_usearch_wide_taxo)
  }
}
################################################################################
