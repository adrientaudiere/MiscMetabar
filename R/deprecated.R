#' Deprecated function(s) in the MiscMetabar package
#'
#' These functions are provided for compatibility with older version of
#' the MiscMetabar package. They may eventually be completely
#' removed.
#' @rdname MiscMetabar-deprecated
#' @name MiscMetabar-deprecated
#' @param ... Parameters to be passed on to the modern version of the function
#' @return Depend on the functions.
#' @export  physeq_graph_test adonis_phyloseq clean_physeq lulu_phyloseq otu_circle biplot_physeq read_phyloseq write_phyloseq sankey_phyloseq summary_plot_phyloseq plot_deseq2_phyloseq plot_edgeR_phyloseq venn_phyloseq ggVenn_phyloseq hill_tuckey_phyloseq hill_phyloseq 
#' @aliases physeq_graph_test adonis_phyloseq clean_physeq lulu_phyloseq otu_circle biplot_physeq read_phyloseq write_phyloseq sankey_phyloseq summary_plot_phyloseq plot_deseq2_phyloseq plot_edgeR_phyloseq venn_phyloseq ggVenn_phyloseq hill_tuckey_phyloseq hill_phyloseq 
#' @section Details:
#' \tabular{rl}{
#'   [graph_test_pq] \tab now a synonym for `physeq_graph_test`\cr
#'   [adonis_pq] \tab now a synonym for `adonis_phyloseq`\cr
#'   [clean_pq] \tab now a synonym for `clean_physeq`\cr
#'   [lulu_pq] \tab now a synonym for `lulu_phyloseq`\cr
#'   [circle_pq] \tab now a synonym for `otu_circle`\cr
#'   [biplot_pq] \tab now a synonym for `biplot_physeq`\cr
#'   [read_pq] \tab now a synonym for `read_phyloseq`\cr
#'   [write_pq] \tab now a synonym for `write_phyloseq`\cr
#'   [sankey_pq] \tab now a synonym for `sankey_phyloseq`\cr
#'   [summary_plot_pq] \tab now a synonym for `summary_plot_phyloseq`\cr
#'   [plot_edgeR_pq] \tab now a synonym for `plot_edgeR_phyloseq`\cr
#'   [plot_deseq2_pq] \tab now a synonym for `plot_deseq2_phyloseq`\cr
#'   [venn_pq] \tab now a synonym for `venn_phyloseq`\cr
#'   [ggvenn_pq] \tab now a synonym for `ggVenn_phyloseq`\cr
#'   [hill_tuckey_pq] \tab now a synonym for `hill_tuckey_phyloseq`\cr
#'   [hill_pq] \tab now a synonym for `hill_phyloseq`\cr
#'   heat_tree_pq \tab now a synonym for `physeq_heat_tree`\cr
#'   [compare_pairs_pq] \tab now a synonym for `multiple_share_bisamples`\cr
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

summary_plot_phyloseq <- function(...) {
  .Deprecated("summary_plot_pq", package = "MiscMetabar")
  summary_plot_pq(...)
}

plot_edgeR_phyloseq <- function(...) {
  .Deprecated("plot_edgeR_pq", package = "MiscMetabar")
  plot_edgeR_pq(...)
}

plot_deseq2_phyloseq <- function(...) {
  .Deprecated("plot_deseq2_pq", package = "MiscMetabar")
  plot_deseq2_pq(...)
}

venn_phyloseq <- function(...) {
  .Deprecated("venn_pq", package = "MiscMetabar")
  venn_pq(...)
}

ggVenn_phyloseq <- function(...) {
  .Deprecated("ggvenn_pq", package = "MiscMetabar")
  ggvenn_pq(...)
}

hill_tuckey_phyloseq <- function(...) {
  .Deprecated("hill_tuckey_pq", package = "MiscMetabar")
  hill_tuckey_pq(...)
}

hill_phyloseq <- function(...) {
  .Deprecated("hill_pq", package = "MiscMetabar")
  hill_pq(...)
}

# physeq_heat_tree <- function(...) {
#   .Deprecated("heat_tree_pq", package = "MiscMetabar")
#   heat_tree_pq(...)
# }

multiple_share_bisamples <- function(...) {
  .Deprecated("compare_pairs_pq", package = "MiscMetabar")
  compare_pairs_pq(...)
}
NULL
