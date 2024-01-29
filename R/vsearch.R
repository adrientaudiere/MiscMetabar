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
#'    vsearch(after filtering taxa based on their number of sequences by the
#'    parameter `min_seq_length` from the [chimera_detection_vs()] function)
#'  - "Select_only_chim" will only select taxa classify as chimera by
#'    vsearch (after filtering taxa based on their number of sequences by the
#'    parameter `min_seq_length` from the [chimera_detection_vs()] function)
#' @param clean_pq (logical; default FALSE) If TRUE, return the phyloseq object
#'   after cleaning using the default parameter of [clean_pq()] function.
#'
#' @param ...
#' @seealso [chimera_detection_vs()]
#' @return
#'
#' - I/ a sequences tables () if object is of class dada, derep, data.frame or
#'   list.
#' - II/ a phyloseq object without (or with if type = 'Select_only_chim')
#'   chimeric taxa
#'
#' @export
#'
#' @examples
#'
#' data_fungi_nochim <- chimera_removal_vs(data_fungi)
#' data_fungi_nochim_16 <- chimera_removal_vs(data_fungi,
#'                                            abskew = 16,
#'                                            min_seq_length = 10)
#' data_fungi_nochim2 <-
#'   chimera_removal_vs(data_fungi, type = "Select_only_non_chim")
#' data_fungi_chimera <-
#'   chimera_removal_vs(data_fungi, type = "Select_only_chim")
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
        chimera_detection_vs(seq2search = colnames(object),
                             nb_seq = colSums(object),
                             ...)
      if (type == "Discard_only_chim") {
        seq_tab_final <-
          object[,!colnames(object) %in% as.character(chim_detect$chimera)]
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
      return(seq_tab_final)

    } else  if (inherits(object, "phyloseq")) {
      chim_detect <-
        chimera_detection_vs(seq2search = refseq(object),
                             nb_seq =  taxa_sums(object),
                             ...)
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
#' @examples
#'   chimera_detection_vs(seq2search = data_fungi@refseq,
#'   nb_seq = taxa_sums(data_fungi))
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
  names(dna_raw) <- paste0("ASV", seq(1, length(seq2search)),
                           ";size=", nb_seq)

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

  Biostrings::writeXStringSet(dna,
                              paste0(tempdir(), "/", "temp.fasta"))

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
      vsearch_arg
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
