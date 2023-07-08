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
    "%>%", ".id", "Ab", "Abundance", "Hill_0", "Hill_1",
    "Hill_2", "X1", "Family", "Genus", "OTU", "Proportion", "Species",
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