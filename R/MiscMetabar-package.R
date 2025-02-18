#' \code{MiscMetabar} package
#'
#' Functions to help analyze and visualize metabarcoding data. Mainly based on
#' the phyloseq and dada2 packages.
#' @name MiscMetabar-package
#' @import ggplot2 phyloseq dada2 dplyr purrr
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".id", "%>%", "Ab", "Abundance", "taxon_names", "beta.div", "calc_taxon_abund",
    "character_method", "Class", "col_tax", "colors", "combn", "complement",
    "e-value", "Family", "Genus", "grid.draw", "grid.layout",
    "group_by", "Guild", "heat_tree", "hill_0", "Hill_0", "hill_1",
    "Hill_1", "hill_2", "Hill_2", "how name", "install_github", "install.packages",
    "LCBD", "log2FoldChange", "logFC", "LVL1", "LVL3", "lwr", "max_Hill",
    "modality", "multcompLetters", "name", "nb_taxa", "nb_seq", "nb_values",
    "ott_id", "OTU", "p.adj", "p.adjust", "plot_layout",
    "plot_layout value", "Proportion", "pushViewport", "Query name",
    "rarefy", "read.delim", "reverse", "reverseComplement", "rgb",
    "rrarefy", "Sample", "Sample_names", "SCBD", "Species", "summarise",
    "tax", "tax_col", "teststat", "Time", "tnrs_match_names", "tol_induced_subtree",
    "update", "upr", "upViewport", "val", "value", "vegdist", "viewport",
    "write.table", "x", "x1", "X1", "x2", "y", "y1", "y2", "ymax",
    "ymin", ".group", "archetype", "nOTUid", "taxon", "total",
    "chim_rm", "condition", "physeq", "seq_tab_Pairs", "nb_samp", "silent",
    "X1_lim1", "X1_lim2", "aicc", "variable", "pos_letters", "alluvium",
    "na_remove", "stratum", "to_lodes_form", "clean_fastq", "clean_sam",
    "samples_names_common", "seq_id", "alpha_hill", "Modality", "Type", "complexity",
    "value_bootstrap", "LearnTaxa", "OrientNucleotides", "RemoveGaps",
    "method", "readDNAStringSet", "Taxa name", "V1", "V2", "pattern_k"
  ))
}

#' @keywords internal
#' @noRd
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats na.exclude na.omit reformulate reorder terms
#' @importFrom lifecycle deprecated
## usethis namespace: end
NULL
