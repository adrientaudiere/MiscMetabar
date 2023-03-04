#' Deprecated function(s) in the MiscMetabar package
#' 
#' These functions are provided for compatibility with older version of
#' the MiscMetabar package.  They may eventually be completely
#' removed.
#' @rdname MiscMetabar-deprecated
#' @name MiscMetabar-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @docType package
#' @export  physeq_graph_test adonis_phyloseq
#' @aliases physeq_graph_test adonis_phyloseq
#' @section Details:
#' \tabular{rl}{
#'   [phyls] \tab now a synonym for [ieua]\cr
#' }
#' 

physeq_graph_test <- function(...) {
  .Deprecated("graph_test_pq",package="yourPackageName")
  graph_test_pq(...)
}

adonis_phyloseq <- function(...) {
  .Deprecated("adonis_pq",package="yourPackageName")
  adonis_pq(...)
}

NULL
