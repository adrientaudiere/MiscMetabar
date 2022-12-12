
################################################################################
#' A wrapper of \code{\link[DECIPHER]{IdTaxa}} which result in a tibble with
#' @description
#' `r lifecycle::badge("experimental")`
#' @param test (required):
#' @param trainingSet (required):
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
