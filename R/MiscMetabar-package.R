#' \code{MiscMetabar} package
#'
#' Functions to help analyze and visualize metabarcoding data. Mainly based on
#' the phyloseq and dada2 packages.
#' @name MiscMetabar-package
#' @docType package
#' @import ggplot2 phyloseq magrittr dada2 dplyr
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".id", "%>%", "Ab", "Abundance", "col_tax", "combn", "complement", "devtools", "e-value", "Family", "Genus", "grid.draw", "grid.layout", "group_by", "Hill_0", "Hill_1", "Hill_2", "install_github", "log2FoldChange", "logFC", "lwr", "max_Hill", "modality", "multcompLetters", "nb_values", "ott_id", "OTU", "Proportion", "pushViewport", "Query name", "rarefy", "reverse", "rgb", "reverseComplement", "rrarefy", "Species", "summarise", "tax", "tax_col", "teststat", "tnrs_match_names", "tol_induced_subtree", "upr", "upViewport", "vegdist", "viewport", "x", "x1", "X1", "x2", "y", "y1", "y2", "ymax", "ymin"
  ))
}

#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats na.exclude na.omit reformulate reorder terms
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
