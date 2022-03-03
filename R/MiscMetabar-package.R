#' \code{MiscMetabar} package
#'
#' Functions to help analyse and visualise metabarcoding data. Mainly based on
#' the phyloseq and dada2 packages.
#' @docType package
#' @name MiscMetabar
#' @import ggplot2 phyloseq magrittr dada2
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "%>%", ".id", "Hill_0", "Hill_1", "Hill_2", "X1",
    "col_tax", "devtools", "grid.draw", "grid.layout",
    "group_by", "install_github", "log2FoldChange",
    "logFC", "lwr", "max_Hill", "modality",
    "multcompLetters", "nb_values", "pushViewport",
    "rarefy", "rrarefy", "summarise", "tax",
    "tax_col", "teststat", "upViewport", "upr",
    "vegdist", "viewport", "x", "x1", "x2", "y",
    "y1", "y2", "ymax", "ymin"
  ))
}

#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
