################################################################################
#' A wrapper of [DECIPHER::IdTaxa()]
#' @description
#' `r lifecycle::badge("experimental")`
#' @param test (required)
#' @param trainingSet (required)
#' @param column_names (optional but often needed): names for the
#' column of the resulted tibble
#' @param ... Additional arguments passed on to \code{\link[DECIPHER]{IdTaxa}}
#'
#' @author Adrien Taudière
#'
#' @return a tibble that can be used for tax_table
#'
#' @export

MM_idtaxa <- function(test, trainingSet, column_names = c("Kingdom_idtaxa", "Phyla_idtaxa", "Class_idtaxa", "Order_idtaxa", "Family_idtaxa", "Genus_idtaxa", "Species_idtaxa", "VT_idtaxa"), ...) {
  idtaxa <- DECIPHER::IdTaxa(test = test, trainingSet = trainingSet, ...)
  col2add <- 8 - lengths(regmatches(idtaxa, gregexpr(";", idtaxa)))
  for (i in seq_along(idtaxa)) {
    idtaxa[i] <-
      paste(
        idtaxa[i],
        paste(as.character(rep(";", each = col2add[i])),
          collapse = ""
        )
      )
  }
  t_idtaxa <- tibble::tibble(data.frame(
    stringr::str_split_fixed(idtaxa, ";", 9)
  ))
  colnames(t_idtaxa) <- column_names
  return(t_idtaxa)
}


#' List fastq files
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @param path path to files (required)
#' @param paired_end do you have paired_end files? (default TRUE)
#' @param pattern a pattern to filter files (passed to list.files function).
#' @param pattern_R1 a pattern to filter R1 files (default "_R1_")
#' @param pattern_R2 a pattern to filter R2 files (default "_R2_")
#' @param nb_files the number of fastq files to list (default FALSE)
#'
#' @return a list of one (single end) or two (paired end) list of files
#'   files are sorted by names (default behavior of `list.files()`)
#' @export
#' @author Adrien Taudière

list_fastq_files <-
  function(path,
           paired_end = TRUE,
           pattern = "fastq",
           pattern_r1 = "_R1_",
           pattern_r2 = "_R2_",
           nb_files = Inf) {
    list_files <- list.files(path, pattern = pattern, full.names = TRUE)
    if (paired_end) {
      fnfs <- sort(list_files[grepl(list_files, pattern = pattern_r1)])
      fnrs <-
        sort(list_files[grepl(list_files, pattern = pattern_r2)])
      if (is.finite(nb_files)) {
        fnfs <- fnfs[1:nb_files]
        fnrs <- fnrs[1:nb_files]
      }
      return(list("fnfs" = fnfs, "fnrs" = fnrs))
    } else {
      fnfs <- sort(list_files[grepl(list_files, pattern = pattern_r1)])
      if (is.finite(nb_files)) {
        fnfs <- fnfs[1:nb_files]
      }
      return(list("fnfs" = fnfs))
    }
  }



#' Rename samples of an otu_table
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @param otu_tab (required) The matrix or otu_table
#' @param names_of_samples (required) The new names of the samples
#' @param taxa_are_rows (default: FALSE) Does the taxa are rows or
#'   columns.
#'
#' @return the matrix with new colnames
#'   (or rownames if `taxa_are_rows` is true)
#'
#' @export
#' @md
#'
#' @author Adrien Taudière
#'
#' @examples
#' data(data_fungi)
#' rename_samples_otu_table(data_fungi@otu_table, as.character(1:185),
#'   taxa_are_rows = T
#' )
rename_samples_otu_table <- function(otu_tab, names_of_samples,
                                     taxa_are_rows = FALSE) {
  if (taxa_are_rows) {
    if (length(names_of_samples) == dim(otu_tab)[1]) {
      rownames(otu_tab) <- names_of_samples
      return(otu_tab)
    } else {
      stop("names_of_samples must have a length equal to the number of samples")
    }
  } else {
    if (length(names_of_samples) == dim(otu_tab)[2]) {
      colnames(otu_tab) <- names_of_samples
      return(otu_tab)
    } else {
      stop("names_of_samples must have a length equal to the number of samples")
    }
  }
}
