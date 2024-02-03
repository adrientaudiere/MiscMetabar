################################################################################
#' List fastq files
#' @description
#' `r lifecycle::badge("maturing")`
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
#' list_fastq_files("extdata")
#' list_fastq_files("extdata", paired_end = FALSE, pattern_R1 = "")
#'
#' @author Adrien Taudière

list_fastq_files <-
  function(path,
           paired_end = TRUE,
           pattern = "fastq",
           pattern_R1 = "_R1_",
           pattern_R2 = "_R2_",
           nb_files = Inf) {
    list_files <- list.files(path, pattern = pattern, full.names = TRUE)
    if (paired_end) {
      fnfs <- sort(list_files[grepl(list_files, pattern = pattern_R1)])
      fnrs <-
        sort(list_files[grepl(list_files, pattern = pattern_R2)])
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
#' `r lifecycle::badge("experimental")`
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
#' data(data_fungi)
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

# ( pattern = "fastq.gz", pattern_R1 = "_R1", pattern_R2 = "_R2")



################################################################################
# WORK IN PROGRESS
#' Remove primers using cutadapt Work In Progress
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @return a list of commande but in The future, will
#'   run the commandes using system
#'
#' @author Adrien Taudière
#'

# cutadapt_remove_primers <- function(
#    path,
#    output_folder = "wo_primers",
#    primer_fw = NULL,
#    primer_rev = NULL,
#    nproc = 1,
#    ...) {
#  cmd <- list()
#
#  if (is.null(primer_rev)) {
#    lff <- list_fastq_files(path, paired_end = FALSE, ...)
#    for (f in lff$fnfs) {
#      cmd[[i]] <- paste0(
#        "source ~/miniconda3/etc/profile.d/conda.sh &&
#    conda activate cutadaptenv && ",
#        "cutadapt --cores=",
#        nproc,
#        " --discard-untrimmed -g '", primer_fw, "' -o ", output_folder,
#        "/", basename(f), " ", f
#      )
#    }
#  } else {
#    lff <- list_fastq_files(path, paired_end = TRUE, ...)
#
#    primer_fw_RC <- dada2:::rc(primer_fw)
#    primer_rev_RC <- dada2:::rc(primer_rev)
#    for (f in lff$fnfs) {
#      cmd[[i]] <- paste0(
#        "source ~/miniconda3/etc/profile.d/conda.sh &&
#    conda activate cutadaptenv && ",
#        "cutadapt -n 2 --cores=",
#        nproc,
#        " --discard-untrimmed -g '", primer_fw,
#        "' -G '", primer_rev,
#        "' -a '", primer_rev_RC,
#        "' -A '", primer_fw_RC,
#        "' -o ", output_folder,
#        "/", basename(f), " -p ", output_folder,
#        "/", gsub("R1", "R2", basename(f)), " ", f, " ",
#        gsub("R1", "R2", f)
#      )
#    }
#  }
#
#  return(cmd)
#  # lapply(cmd, system)
# }
################################################################################


################################################################################
#' A wrapper of the function [dada2::filterAndTrim()] to use in
#'   [targets](https://books.ropensci.org/targets/) pipeline
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' This function filter and trim (with parameters passed on to
#'   [dada2::filterAndTrim()] function) forward sequences or paired end
#'   sequence if 'rev' parameter is set. It return the list of files to
#'   subsequent analysis in a targets pipeline.
#'
#'
#' @param fw (required) a list of forward fastq files
#' @param rev a list of reverse fastq files for paired end trimming
#' @param output_fw Path to output folder for forward files. By default,
#'   this function will create
#'   a folder "output/filterAndTrim_fwd" in the current working directory.
#' @param output_rev Path to output folder for reverse files. By default,
#'   this function will create
#'   a folder "output/filterAndTrim_fwd" in the current working directory.
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
#'
#' @seealso [dada2::filterAndTrim()]
filter_trim <-
  function(fw = NULL,
           rev = NULL,
           output_fw = paste(getwd(), "/output/filterAndTrim_fwd", sep = ""),
           output_rev = paste(getwd(), "/output/filterAndTrim_rev", sep = ""),
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
        dada2::filterAndTrim(
          filt = output_fw,
          fwd = fw, ...
        )
        return(output_fw)
      } else {
        dada2::filterAndTrim(
          filt = output_fw,
          filt.rev = output_rev,
          fwd = fw,
          rev = rev,
          ...
        )
        return(list("fw" = output_fw, "rv" = output_rev))
      }
    }
  }
################################################################################


################################################################################
#' Load sample data from file and rename samples using names of samples and an
#'   optional order
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @param file_path (required) a path to the sample_data file
#' @param names_of_samples (required) a vector of sample names
#' @param samples_order Optional numeric vector to sort sample names
#' @param ... Other arguments passed on to [utils::read.delim()] function.
#'
#' @return A data.frame from file_path and new names
#' @export
#'
#' @author Adrien Taudière
#'
#' @seealso [rename_samples()]
sample_data_with_new_names <- function(file,
                                       names_of_samples,
                                       samples_order,
                                       ...) {
  # For use with old syntax using condition in place of samples_order
  if (is.null(samples_order) && !missing(condition)) {
    samples_order <- condition
  }
  samdata_interm <- read.delim(file, ...)
  samdata_renamed <- rename_samples(samdata_interm, names_of_samples)
  if (missing(samples_order)) {
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
#' `r lifecycle::badge("maturing")`
#'
#' @param phyloseq_component (required) one of otu_table or sam_data slot of a
#'   phyloseq-class object
#' @param names_of_samples (required) A vector of samples names
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
                           names_of_samples) {
  if (length(names_of_samples) != length(sample_names(phyloseq_component))) {
    stop("Names_of_samples must have a length equal to the number of samples.")
  }

  new_pq_component <- phyloseq_component
  sample_names(new_pq_component) <- names_of_samples
  return(new_pq_component)
}
################################################################################
