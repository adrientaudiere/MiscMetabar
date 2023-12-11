#' Fungal OTU in phyloseq format
#'
#' @format A physeq object containing 1420 taxa with references sequences
#' described by 14 taxonomic ranks and 185 samples described by 7 sample variables:
#' - *X*: the name of the fastq-file
#' - *Sample_names*: the names of ... the samples
#' - *Treename*: the name of an tree
#' - *Sample_id*: identifier for each sample
#' - *Height*: height of the sample in the tree
#' - *Diameter*: diameter of the trunk
#' - *Time*: time since the dead of the tree
"data_fungi"

#' Fungal OTU in phyloseq format
#'
#' It is a subset of the data_fungi dataset including only taxa with information
#' at the species level
#'
#' Obtain using `data_fungi_sp_known <- subset_taxa(data_fungi,
#'   !is.na(data_fungi@tax_table[,"Species"]))`
#'
#' @format A physeq object containing 651 taxa with references sequences
#' described by 14 taxonomic ranks and 185 samples described by 7 sample variables:
#' - *X*: the name of the fastq-file
#' - *Sample_names*: the names of ... the samples
#' - *Treename*: the name of an tree
#' - *Sample_id*: identifier for each sample
#' - *Height*: height of the sample in the tree
#' - *Diameter*: diameter of the trunk
#' - *Time*: time since the dead of the tree
"data_fungi_sp_known"

#' This tutorial explore the dataset from Tengeler et al. (2020) available in the `mia` package.
#' obtain using `mia::makePhyloseqFromTreeSE(Tengeler2020)` 
#' 
#' This is a phyloseq version of the Tengeler2020 dataset.
#' 
#' Tengeler2020 includes gut microbiota profiles of 27 persons with ADHD.
#' A standard bioinformatic and statistical analysis done to demonstrate that altered microbial 
#' composition could be a driver of altered brain structure and function and concomitant changes 
#' in the animals behavior. This was investigated by colonizing young, male, 
#' germ-free C57BL/6JOlaHsd mice with microbiota from individuals with and without ADHD.
#'
#' @format 
#' A phyloseq object
#' 
"Tengeler2020_pq"