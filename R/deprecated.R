#' Deprecated function(s) in the MiscMetabar package
#'
#' These functions are provided for compatibility with older version of
#' the MiscMetabar package. They may eventually be completely
#' removed.
#' @rdname MiscMetabar-deprecated
#' @name MiscMetabar-deprecated
#' @param ... Parameters to be passed to the modern version of the function
#' @docType package
#' @export  physeq_graph_test adonis_phyloseq clean_physeq lulu_phyloseq otu_circle biplot_physeq read_phyloseq write_phyloseq sankey_phyloseq
#' @aliases physeq_graph_test adonis_phyloseq clean_physeq lulu_phyloseq otu_circle biplot_physeq read_phyloseq write_phyloseq sankey_phyloseq
#' @section Details:
#' \tabular{rl}{
#'   [graph_test_pq] \tab now a synonym for [physeq_graph_test]\cr
#'   [adonis_pq] \tab now a synonym for [adonis_phyloseq]\cr
#'   [clean_pq] \tab now a synonym for [clean_physeq]\cr
#'   [lulu_pq] \tab now a synonym for [lulu_phyloseq]\cr
#'   [circle_pq] \tab now a synonym for [otu_circle]\cr
#'   [biplot_pq] \tab now a synonym for [biplot_physeq]\cr
#'   [read_pq] \tab now a synonym for [read_phyloseq]\cr
#'   [write_pq] \tab now a synonym for [write_phyloseq]\cr
#'   [sankey_pq] \tab now a synonym for [sankey_phyloseq]\cr
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

read_phyloseq <- function(...) {
  .Deprecated("read_pq", package = "MiscMetabar")
  read_pq(...)
}

write_phyloseq <- function(...) {
  .Deprecated("write_pq", package = "MiscMetabar")
  write_pq(...)
}

sankey_phyloseq <- function(...) {
  .Deprecated("sankey_pq", package = "MiscMetabar")
  sankey_pq(...)
}

NULL