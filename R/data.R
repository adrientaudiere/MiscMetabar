#' Some fungal OTU in phyloseq format
#'
#' @format A physeq object containing 1420 taxa with references sequences
#' described by 14 taxonomic ranks
#' and 185 samples described by 7 sample variables:
"data_fungi"

#' Some fungal OTU in phyloseq format with known species name
#'
#' A subset of the data_fungi dataset including only taxa with information
#' at the species level
#' 
#' Obtain using `data_fungi_sp_known <- subset_taxa(data_fungi, !is.na(data_fungi@tax_table[,"Species"]))`
#'
#' @format A physeq object containing 651 taxa with references sequences
#' described by 14 taxonomic ranks
#' and 185 samples described by 7 sample variables:
"data_fungi_sp_known"
