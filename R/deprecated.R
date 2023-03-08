#' Deprecated function(s) in the MiscMetabar package
#'
#' These functions are provided for compatibility with older version of
#' the MiscMetabar package. They may eventually be completely
#' removed.
#' @rdname MiscMetabar-deprecated
#' @name MiscMetabar-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @docType package
#' @export  physeq_graph_test adonis_phyloseq clean_physeq lulu_phyloseq otu_circle biplot_physeq
#' @aliases physeq_graph_test adonis_phyloseq clean_physeq lulu_phyloseq otu_circle biplot_physeq
#' @section Details:
#' \tabular{rl}{
#'   [graph_test_pq] \tab now a synonym for [physeq_graph_test]\cr
#'   [adonis_pq] \tab now a synonym for [adonis_phyloseq]\cr
#'   [clean_pq] \tab now a synonym for [clean_physeq]\cr
#'   [lulu_pq] \tab now a synonym for [lulu_phyloseq]\cr
#'   [circle_pq] \tab now a synonym for [otu_circle]\cr
#'   [biplot_pq] \tab now a synonym for [biplot_physeq]\cr
#' }
#'

physeq_graph_test <- function(...) {
  .Deprecated("graph_test_pq", package = "MiscMetabar")
  graph_test_pq(...)
}

adonis_phyloseq <- function(...) {
  .Deprecated("adonis_pq", package = "MiscMetabar")
  adonis_pq(...)
}

clean_physeq <- function(...) {
  .Deprecated("clean_pq", package = "MiscMetabar")
  clean_pq(...)
}


lulu_phyloseq <- function(...) {
  .Deprecated("lulu_pq", package = "MiscMetabar")
  lulu_pq(...)
}

otu_circle <- function(...) {
  .Deprecated("circle_pq", package = "MiscMetabar")
  circle_pq(...)
}

biplot_physeq <- function(...) {
  .Deprecated("biplot_pq", package = "MiscMetabar")
  biplot_pq(...)
}

NULL