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
#' @usage data(data_fungi)
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
#' @usage data(data_fungi_sp_known)
"data_fungi_sp_known"


#' Fungal OTU in phyloseq format
#'
#' It is a subset of the data_fungi dataset including only Basidiomycota
#'   with more than 5000 sequences.
#'
#' Obtain using `data_fungi_mini <- subset_taxa(data_fungi, Phylum == "Basidiomycota")`
#' and then `data_fungi_mini <-   subset_taxa_pq(data_fungi_mini, colSums(data_fungi_mini@otu_table) > 5000)`
#'
#' @format A physeq object containing 45 taxa with references sequences
#' described by 14 taxonomic ranks and 137 samples described by 7 sample variables:
#' - *X*: the name of the fastq-file
#' - *Sample_names*: the names of ... the samples
#' - *Treename*: the name of an tree
#' - *Sample_id*: identifier for each sample
#' - *Height*: height of the sample in the tree
#' - *Diameter*: diameter of the trunk
#' - *Time*: time since the dead of the tree
#' @usage data(data_fungi_mini)
#' @format A physeq object containing 45 taxa with references sequences
#' described by 14 taxonomic ranks and 137 samples described by 7 sample variables:
#' - *X*: the name of the fastq-file
#' - *Sample_names*: the names of ... the samples
#' - *Treename*: the name of an tree
#' - *Sample_id*: identifier for each sample
#' - *Height*: height of the sample in the tree
#' - *Diameter*: diameter of the trunk
#' - *Time*: time since the dead of the tree
#' @usage data(data_fungi_mini)
"data_fungi_mini"


#' This tutorial explores the dataset from Tengeler et al. (2020) available in the `mia` package.
#' obtained using `mia::makePhyloseqFromTreeSE(Tengeler2020)`
#'
#' This is a phyloseq version of the Tengeler2020 dataset.
#'
#' Tengeler2020 includes gut microbiota profiles of 27 persons with ADHD.
#' A standard bioinformatic and statistical analysis done to demonstrate that altered microbial
#' composition could be a driver of altered brain structure and function and concomitant changes
#' in the animals behavior. This was investigated by colonizing young, male,
#' germ-free C57BL/6JOlaHsd mice with microbiota from individuals with and without ADHD.
#'
#' Tengeler, A.C., Dam, S.A., Wiesmann, M. et al. Gut microbiota from persons with
#' attention-deficit/hyperactivity disorder affects the brain in mice. Microbiome 8, 44 (2020).
#' https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00816-x
#' @format
#' A phyloseq object
#' @usage data(Tengeler2020_pq)
"Tengeler2020_pq"


#' Default patterns for unwanted taxonomic values
#'
#' @description
#' A named character vector of regular expressions used to identify common
#' problematic values in taxonomy tables. Each element is a regex pattern;
#' names provide human-readable descriptions.
#'
#' Used as the default `replace_to_NA` argument in [verify_tax_table()] and
#' can be reused by other pqverse packages (e.g. `dbpq::count_unwanted_tax()`).
#'
#' @format A named character vector with 17 elements:
#' \describe{
#'   \item{NA-like (NA, NaN, nan)}{`"^[Nn][Aa][Nn]?$"`}
#'   \item{NA-like (N/A, n/a)}{`"^[Nn]/[Aa]$"`}
#'   \item{None / none}{`"^[Nn]one$"`}
#'   \item{empty string}{`"^$"`}
#'   \item{whitespace only}{`"^\\\\s+$"`}
#'   \item{unclassified}{`"[Uu]nclassified"`}
#'   \item{unknown}{`"[Uu]nknown"`}
#'   \item{unidentified}{`"[Uu]nidentified"`}
#'   \item{uncultured}{`"[Uu]ncultured"`}
#'   \item{incertae sedis}{`"[Ii]ncertae[_\\\\s]?[Ss]edis"`}
#'   \item{metagenome}{`"^[Mm]etagenome$"`}
#'   \item{environmental}{`"^[Ee]nvironmental"`}
#'   \item{empty QIIME-style rank}{`"^[kpcofgs]__$"`}
#'   \item{unknown species (_sp prefix)}{`"^_sp"`}
#'   \item{unknown species (_species prefix)}{`"^_species"`}
#'   \item{unknown cluster (MMseqs2)}{`"_uc$"`}
#'   \item{unknown ranks (PR2 database)}{`"__X+$"`}
#' }
#'
#' @export
#' @seealso [verify_tax_table()]
#' @examples
#' unwanted_tax_patterns
#' # Use with grepl to check a value
#' any(vapply(
#'   unwanted_tax_patterns,
#'   \(pat) grepl(pat, "unclassified"),
#'   logical(1)
#' ))
unwanted_tax_patterns <- c(
  "NA-like (NA, NaN, nan)" = "^[Nn][Aa][Nn]?$",
  "NA-like (N/A, n/a)" = "^[Nn]/[Aa]$",
  "None / none" = "^[Nn]one$",
  "empty string" = "^$",
  "whitespace only" = "^\\s+$",
  "unclassified" = "[Uu]nclassified",
  "unknown" = "[Uu]nknown",
  "unidentified" = "[Uu]nidentified",
  "uncultured" = "[Uu]ncultured",
  "incertae sedis" = "[Ii]ncertae[_\\s]?[Ss]edis",
  "metagenome" = "^[Mm]etagenome$",
  "environmental" = "^[Ee]nvironmental",
  "empty QIIME-style rank" = "^[kpcofgs]__$",
  "unknown species (_sp prefix)" = "^_sp",
  "unknown species (_species prefix)" = "^_species",
  "unknown cluster (_uc prefix, e.g. MMseqs2 assignation)" = "_uc$",
  "unknown ranks (_X, _XX, ... prefix e.g. PR2 database)" = "__X+$"
)
