################################################################################
#' List fastq files
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Useful for targets bioinformatic pipeline.
#'
#' @param path path to files (required)
#' @param paired_end do you have paired_end files? (default TRUE)
#' @param pattern a pattern to filter files (passed on to list.files function).
#' @param pattern_R1 a pattern to filter R1 files (default "_R1_")
#' @param pattern_R2 a pattern to filter R2 files (default "_R2_")
#' @param nb_files the number of fastq files to list (default FALSE)
#'
#' @return a list of one (single end) or two (paired end) list of files
#'   files are sorted by names (default behavior of `list.files()`)
#' @export
#'
#' @examples
#' list_fastq_files(system.file("extdata", package = "MiscMetabar"))
#' list_fastq_files(system.file("extdata", package = "MiscMetabar"),
#'   paired_end = FALSE, pattern_R1 = ""
#' )
#'
#' @author Adrien Taudière

list_fastq_files <-
  function(path,
           paired_end = TRUE,
           pattern = "fastq",
           pattern_R1 = "_R1_",
           pattern_R2 = "_R2_",
           nb_files = Inf) {
    if (length(list.files(path)) == 0) {
      stop("There is no files in the folder ", path)
    }
    list_files <- list.files(path, pattern = pattern, full.names = TRUE)
    if (length(list_files) == 0) {
      stop("None file in the folder ", path, " match the pattern ", pattern)
    }
    if (paired_end) {
      fnfs <- sort(list_files[grepl(list_files, pattern = pattern_R1)])
      if (length(fnfs) == 0) {
        stop("None file in the folder ", path, " match the pattern_R1 ", pattern_R1)
      }
      fnrs <-
        sort(list_files[grepl(list_files, pattern = pattern_R2)])
      if (length(fnrs) == 0) {
        stop("None file in the folder ", path, " match the pattern_R2 ", pattern_R2)
      }
      if (is.finite(nb_files)) {
        fnfs <- fnfs[seq(1, nb_files)]
        fnrs <- fnrs[seq(1, nb_files)]
      }
      return(list("fnfs" = fnfs, "fnrs" = fnrs))
    } else {
      fnfs <- sort(list_files[grepl(list_files, pattern = pattern_R1)])
      if (is.finite(nb_files)) {
        fnfs <- fnfs[seq(1, nb_files)]
      }
      return(list("fnfs" = fnfs))
    }
  }
################################################################################


################################################################################
#' Rename samples of an otu_table
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful for targets bioinformatic pipeline.
#'
#' @inheritParams clean_pq
#' @param names_of_samples (required) The new names of the samples
#'
#' @return the matrix with new colnames (or rownames if `taxa_are_rows` is true)
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#'
#' rename_samples_otu_table(data_fungi, as.character(seq_along(sample_names(data_fungi))))
#'
rename_samples_otu_table <- function(physeq, names_of_samples) {
  otu_tab <- physeq@otu_table
  tax_in_row <- taxa_are_rows(physeq)
  if (tax_in_row) {
    if (length(names_of_samples) == dim(otu_tab)[2]) {
      colnames(otu_tab) <- names_of_samples
      return(otu_tab)
    } else {
      stop("names_of_samples must have a length equal to the number of samples")
    }
  } else {
    if (length(names_of_samples) == dim(otu_tab)[1]) {
      rownames(otu_tab) <- names_of_samples
      return(otu_tab)
    } else {
      stop("names_of_samples must have a length equal to the number of samples")
    }
  }
}
################################################################################

################################################################################
#' A wrapper of the function [dada2::filterAndTrim()] to use in
#'   [targets](https://books.ropensci.org/targets/) pipeline
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' This function filter and trim (with parameters passed on to
#'   [dada2::filterAndTrim()] function) forward sequences or paired end
#'   sequence if 'rev' parameter is set. It return the list of files to
#'   subsequent analysis in a targets pipeline.
#'
#' @param fw (required) a list of forward fastq files
#' @param rev a list of reverse fastq files for paired end trimming
#' @param output_fw Path to output folder for forward files. By default,
#'   this function will create
#'   a folder "output/filterAndTrim_fwd" in the current working directory.
#' @param output_rev Path to output folder for reverse files. By default,
#'   this function will create
#'   a folder "output/filterAndTrim_fwd" in the current working directory.
#' @param return_a_vector (logical, default FALSE) If true, the return is
#'   a vector of path (usefull when used with
#'   targets::tar_targets(..., format="file"))
#' @param ... Other parameters passed on to [dada2::filterAndTrim()] function.
#'
#' @return A list of files. If rev is set, will return a list of two lists.
#'  The first list is a list of forward files, and the second one
#'  is a list of reverse files.
#' @export
#'
#' @examples
#' testFastqs_fw <- c(
#'   system.file("extdata", "sam1F.fastq.gz", package = "dada2"),
#'   system.file("extdata", "sam2F.fastq.gz", package = "dada2")
#' )
#' testFastqs_rev <- c(
#'   system.file("extdata", "sam1R.fastq.gz", package = "dada2"),
#'   system.file("extdata", "sam2R.fastq.gz", package = "dada2")
#' )
#'
#' filt_fastq_fw <- filter_trim(testFastqs_fw, output_fw = tempdir())
#' derep_fw <- derepFastq(filt_fastq_fw[1])
#' derep_fw
#'
#' filt_fastq_pe <- filter_trim(testFastqs_fw,
#'   testFastqs_rev,
#'   output_fw = tempdir("fw"),
#'   output_rev = tempdir("rev")
#' )
#' derep_fw_pe <- derepFastq(filt_fastq_pe[[1]])
#' derep_rv_pe <- derepFastq(filt_fastq_pe[[2]])
#' derep_fw_pe
#' derep_rv_pe
#' @author Adrien Taudière
#'
#' @seealso [dada2::filterAndTrim()]
filter_trim <-
  function(fw = NULL,
           rev = NULL,
           output_fw = file.path(paste(getwd(), "/output/filterAndTrim_fwd", sep = "")),
           output_rev = file.path(paste(getwd(), "/output/filterAndTrim_rev", sep = "")),
           return_a_vector = FALSE,
           ...) {
    if (length(fw) == 1) {
      # This case with one file is to create a folder instead of only one file

      if (!is.null(rev)) {
        dada2::filterAndTrim(
          filt = paste0(output_fw, "interm"),
          filt.rev = paste0(output_rev, "interm"),
          fwd = fw,
          rev = rev,
          ...
        )
        dir.create(output_rev)
        dir.create(output_fw)
        file.rename(
          paste0(output_rev, "interm"),
          paste0(output_rev, "/", basename(rev))
        )
        file.rename(
          paste0(output_fw, "interm"),
          paste0(output_fw, "/", basename(fw))
        )

        return(list("fw" = output_fw, "rv" = output_rev))
      } else {
        dada2::filterAndTrim(
          filt = paste0(output_fw, "interm"),
          fwd = fw,
          ...
        )
        dir.create(output_fw)

        file.rename(
          paste0(output_fw, "interm"),
          paste0(output_fw, "/", basename(fw))
        )
        return(output_fw)
      }
    } else {
      if (is.null(rev)) {
        dada2::filterAndTrim(filt = output_fw, fwd = fw, ...)
        return(output_fw)
      } else {
        dada2::filterAndTrim(
          filt = output_fw,
          filt.rev = output_rev,
          fwd = fw,
          rev = rev,
          ...
        )
        if (return_a_vector) {
          return(c("fw" = output_fw, "rv" = output_rev))
        } else {
          return(list("fw" = output_fw, "rv" = output_rev))
        }
      }
    }
  }
################################################################################


################################################################################
#' Load sample data from file and rename samples using names of samples and an
#'   optional order
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Useful for targets bioinformatic pipeline.
#'
#' @param file_path (required) a path to the sample_data file
#' @param names_of_samples (required) a vector of sample names
#' @param samples_order Optional numeric vector to sort sample names
#' @param ... Additional arguments passed on to [utils::read.delim()] function.
#'
#' @return A data.frame from file_path and new names
#' @export
#'
#' @examples
#' sam_file <- system.file("extdata", "sam_data.csv", package = "MiscMetabar")
#' sample_data_with_new_names(sam_file, paste0("Samples_", seq(1, 185)))
#'
#' @author Adrien Taudière
#'
#' @seealso [rename_samples()]
sample_data_with_new_names <- function(file_path,
                                       names_of_samples,
                                       samples_order = NULL,
                                       ...) {
  samdata_interm <- sample_data(read.delim(file_path, ...))
  samdata_renamed <- rename_samples(samdata_interm, names_of_samples)
  if (is.null(samples_order)) {
    samdata <- samdata_renamed
  } else {
    samdata <- samdata_renamed[samples_order, ]
  }
  return(samdata)
}
################################################################################

################################################################################
#' Rename the samples of a phyloseq slot
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Useful for targets bioinformatic pipeline.
#'
#' @param phyloseq_component (required) one of otu_table or sam_data slot of a
#'   phyloseq-class object
#' @param names_of_samples (required) A vector of samples names
#' @param taxa_are_rows (default to FALSE) see ?phyloseq for details
#'
#' @return The otu_table or the sam_data slot with new samples names
#' @export
#' @author Adrien Taudière
#' @examples
#' otutab <- rename_samples(
#'   data_fungi@otu_table,
#'   paste0("data_f", sample_names(data_fungi))
#' )
#' otutab2 <- rename_samples(
#'   clean_pq(data_fungi,
#'     force_taxa_as_rows = TRUE
#'   )@otu_table,
#'   paste0("data_f", sample_names(data_fungi))
#' )
#' samda <- rename_samples(
#'   data_fungi@sam_data,
#'   paste0("data_f", sample_names(data_fungi))
#' )
rename_samples <- function(phyloseq_component,
                           names_of_samples,
                           taxa_are_rows = FALSE) {
  if (is.null(sample_names(phyloseq_component)) &&
    inherits(phyloseq_component, "matrix")) {
    phyloseq_component <- otu_table(phyloseq_component, taxa_are_rows = taxa_are_rows)
  }
  if (length(names_of_samples) != length(sample_names(phyloseq_component))) {
    stop("Names_of_samples must have a length equal to the number of samples.")
  }

  new_pq_component <- phyloseq_component
  sample_names(new_pq_component) <- names_of_samples
  return(new_pq_component)
}
################################################################################


################################################################################
#' Match sample names from sam_data and fastq files
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful for targets bioinformatic pipeline.
#'
#' @param path_sam_data (Required) Path to sample data file.
#' @param sample_col_name (Required) The name of the column defining sample names in the sample data file.
#' @param path_raw_seq (Required) Path to the folder containing fastq files
#' @param pattern_remove_sam_data If not null, describe the pattern that will be deleted from
#'   sam_data samples names.
#' @param pattern_remove_fastq_files If not null, describe the pattern that will be deleted from
#'   fastq files names.
#' @param verbose (logical, default TRUE) If TRUE, print some additional messages.
#' @param remove_undocumented_fastq_files (logical, default FALSE) If set to TRUE
#'   fastq files not present in sam_data are removed from your folder.
#'   Keep a copy of those files somewhere before.
#' @param prefix Add a prefix to new samples names (ex. prefix = "samp")
#' @param ... Other parameters passed on to [utils::read.csv()] function.
#' @return A list of two objects :
#'   - $sam_names_matching is a tibble of corresponding samples names
#'   - $sam_data is a sample data files including only matching sample names
#' @importFrom utils read.csv
#' @export
#' @author Adrien Taudière
sam_data_matching_names <- function(path_sam_data,
                                    sample_col_name,
                                    path_raw_seq,
                                    pattern_remove_sam_data = NULL,
                                    pattern_remove_fastq_files = NULL,
                                    verbose = TRUE,
                                    remove_undocumented_fastq_files = FALSE,
                                    prefix = NULL,
                                    ...) {
  sam_d <- read.csv(path_sam_data, ...)
  names_sam_data <- sam_d[[sample_col_name]]
  names_fastq_files_fullpath <- list.files(path_raw_seq)
  names_fastq_files <- basename(names_fastq_files_fullpath)

  if (!is.null(pattern_remove_sam_data)) {
    names_sam_data_clean <- gsub(pattern_remove_sam_data, "", names_sam_data)
    if (sum(duplicated(names_sam_data_clean)) > 0) {
      stop(
        "Their is duplicates in samples names from sam_data when removing pattern_remove_sam_data",
        names_sam_data_clean[duplicated(names_sam_data_clean)]
      )
    }
  } else {
    names_sam_data_clean <- names_sam_data
  }

  if (!is.null(pattern_remove_fastq_files)) {
    names_fastq_files_clean <- gsub(pattern_remove_fastq_files, "", names_fastq_files)
  } else {
    names_fastq_files_clean <- names_fastq_files
  }

  tib_sam_data <- tibble(
    clean_sam = as.character(names_sam_data_clean),
    raw_sam = as.character(names_sam_data)
  )

  tib_fastq <- tibble(
    clean_fastq = as.character(names_fastq_files_clean),
    raw_fastq = as.character(names_fastq_files),
    raw_fastq_full_path = names_fastq_files_fullpath
  )

  tib_j <- full_join(tib_fastq, tib_sam_data, by = join_by(clean_fastq == clean_sam)) |>
    rename(common_names = clean_fastq)

  if (sum(is.na(tib_j$raw_fastq)) > 0) {
    message(
      sum(is.na(tib_j$raw_fastq)),
      " samples in sam_data files are not present in fastq_files"
    )
    if (verbose) {
      warning(
        paste(tib_j$raw_sam[is.na(tib_j$raw_fastq)], collapse = " "),
        " not matching names from sam_data file."
      )
    }
  }

  if (sum(is.na(tib_j$raw_sam)) > 0) {
    message(
      sum(is.na(tib_j$raw_sam)),
      " samples in fastq files are not present in sam_data"
    )
    if (verbose) {
      warning(
        paste(tib_j$raw_fastq[is.na(tib_j$raw_sam)], collapse = " "),
        " not matching names from fastq files."
      )
    }
    if (remove_undocumented_fastq_files) {
      if (verbose) {
        warning(
          "Files ",
          tib_j$raw_fastq_full_path[is.na(tib_j$raw_sam)],
          "will be removed from folder",
          path_raw_seq
        )
      }
      unlink(tib_j$raw_fastq_full_path[is.na(tib_j$raw_sam)])
    }
  }

  if (!is.null(prefix)) {
    tib_j$common_names <- paste0(prefix, tib_j$common_names)
  }

  if (is.null(pattern_remove_sam_data)) {
    sam_d_new <- sam_d |>
      dplyr::filter(.data[[sample_col_name]] %in% tib_j$raw_sam[!is.na(tib_j$raw_fastq)]) |>
      dplyr::mutate("samples_names_common" = paste0(prefix, .data[[sample_col_name]])) |>
      relocate(samples_names_common)
  } else {
    sam_d_new <- sam_d |>
      dplyr::filter(.data[[sample_col_name]] %in% tib_j$raw_sam[!is.na(tib_j$raw_fastq)]) |>
      dplyr::mutate("samples_names_common" = paste0(prefix, gsub(pattern_remove_sam_data, "", .data[[sample_col_name]]))) |>
      relocate(samples_names_common)
  }
  return(list("sam_names_matching" = tib_j, "sam_data" = sam_d_new))
}
################################################################################
