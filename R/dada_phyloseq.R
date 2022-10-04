if (getRversion() >= "2.15.1") utils::globalVariables(c("."))

################################################################################
#' Add dna in `refseq` slot of a physeq object using taxa names and renames taxa
#'   using ASV_1, ASV_2, …
#'
#' `r lifecycle::badge("stable")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object obtained
#'   using the `dada2` package
#'
#' @return A new \code{\link{phyloseq-class}} object with `refseq` slot and new
#'   taxa names
#' @export

add_dna_to_phyloseq <- function(physeq) {
  if (!methods::validObject(physeq) || class(physeq) != "phyloseq") {
    stop("The physeq argument is not a valid phyloseq object.")
  } else {
    dna <- Biostrings::DNAStringSet(phyloseq::taxa_names(physeq))
    names(dna) <- phyloseq::taxa_names(physeq)
    physeq <- phyloseq::merge_phyloseq(physeq, dna)
    phyloseq::taxa_names(physeq) <-
      paste0("ASV_", seq(phyloseq::ntaxa(physeq)))
  }
  return(physeq)
}
################################################################################



################################################################################
#'  Clean phyloseq object by removing empty samples and taxa
#'  In addition, this function check for discrepancy (and rename) between
#' (i) taxa names in refseq, taxonomy table and otu_table and between
#' (ii) sample names in sam_data and otu_table.
#'
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object obtained
#'   using the `dada2` package
#'
#' @param remove_empty_samples (logical, default TRUE) Do you want to remove samples without sequences (this is done after removing empty taxa)
#' @param remove_empty_taxa (logical, default TRUE) Do you want to remove taxa without sequences (this is done before removing empty samples)
#' @param clean_samples_names
#' @return A new \code{\link{phyloseq-class}} object
#' @export
clean_physeq <-  function(physeq,
                          remove_empty_samples = TRUE,
                          remove_empty_taxa = TRUE,
                          clean_samples_names =TRUE
                          ){
  if (clean_samples_names) {
    if (!is.null(physeq@refseq)){
      if (sum(!names(physeq@refseq) %in% taxa_names(physeq)) > 0){
        names(physeq@refseq) <- taxa_names(physeq)
        message("Change the samples names in refseq slot")
      }
    }
    if (!is.null(physeq@tax_table)){
      if (sum(!rownames(physeq@tax_table) %in% taxa_names(physeq)) > 0){
        rownames(physeq@tax_table) <- taxa_names(physeq)
        message("Change the taxa names in tax_table slot")
      }
    }

    if (!is.null(physeq@sam_data)){
      if (sum(!rownames(physeq@sam_data) %in% sample_names(physeq)) > 0){
        rownames(physeq@sam_data) <- sample_names(physeq)
        message("Change the samples names in sam_data slot")
      }
    }
  }

  if (sum(grepl("^0", "", sample_names(physeq)) )){
    message("At least one sample name start with a zero. That can be a problem for some phyloseq functions such as plot_bar and psmelt.")
  }

  new_physeq <- physeq

  if (remove_empty_taxa) {
    if (sum(taxa_sums(new_physeq) == 0) > 0) {
    # new_otu_table <- otu_table(new_physeq, taxa_are_rows =T)[,taxa_sums(new_physeq) > 0]
    # new_tax_table <- tax_table(new_physeq)[taxa_sums(new_physeq) > 0,]
    # new_physeq <- merge_phyloseq(new_otu_table, new_tax_table, physeq)
      new_physeq <- subset_taxa(physeq, taxa_sums(physeq)>0)
    }
  }
  if (remove_empty_samples) {
    if (sum(sample_sums(new_physeq) == 0) > 0) {
      new_physeq <- subset_samples(new_physeq, sample_sums(physeq) > 0)
    }
  }
  message(paste("Supress", ntaxa(physeq)-ntaxa(new_physeq), "taxa and",
          nsamples(physeq)-nsamples(new_physeq),
          "samples.")
         )
  return(new_physeq)
}






################################################################################
#' Track the number of reads (= sequences), samples and cluster (e.g. ASV)
#' from various objects including dada-class and derep-class.
#'
#' `r lifecycle::badge("maturing")`
#'
#'  * List of fastq and fastg.gz files -> nb of reads and samples
#'  * List of dada-class -> nb of reads, clusters (ASV) and samples
#'  * List of derep-class -> nb of reads, clusters (unique sequences)
#'    and samples
#'  * Matrix of samples x clusters (e.g. otu_table) -> nb of reads,
#'    clusters and samples
#'  * Phyloseq-class -> nb of reads, clusters and samples
#'
#' @param list_of_objects (required): A list of objects
#' @param obj_names (default: NULL) :
#'   A list of names corresponding to the list of objects
#' @param clean_physeq (Logical, default = FALSE) If true, empty samples and empty ASV are discarded before clustering.
#'
#' @return The number of sequences, clusters (e.g. OTUs, ASVs) and samples for
#'   each object.
#' @export

track_wkflow <- function(list_of_objects, obj_names = NULL, clean_physeq = FALSE) {
  message("Compute the number of sequences")
  if (!is.null(obj_names)) {
    names(list_of_objects) <- obj_names
  }

  if (clean_physeq) {
    for (i in seq_along(list_of_objects)) {
      if (inherits(list_of_objects[[i]], "phyloseq")) {
        list_of_objects[[i]] <- clean_physeq(list_of_objects[[i]])
      }
    }
  }

  track_nb_seq_per_obj <-
    pbapply::pblapply(list_of_objects, function(object) {
      message(paste("Start object of class:", class(object), sep = " "))
      if (inherits(object, "phyloseq")) {
        sum(object@otu_table)
      } else if (inherits(object, "matrix")) {
        sum(object, na.rm = TRUE)
      } else if (is.character(object[1]) &&
        length(object[1]) == 1 &&
        file.exists(object[1])) {
        if (summary(file(object[[1]]))$class == "gzfile") {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("zcat ", x, " | grep -c '^+$'", sep = ""),
              intern = TRUE
            ))
          })
        } else if (grepl("\\.fastq$", object[1])) {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("cat ", x, " | grep -c '^+$'", sep = ""),
              intern = TRUE
            ))
          })
        } else {
          stop("Files must be either gzfile or .fastq")
        }
      } else {
        pbapply::pbsapply(object, function(x) {
          sum(dada2::getUniques(x, silence = TRUE))
        })
      }
    })
  track_nb_seq_per_obj <-
    pbapply::pblapply(track_nb_seq_per_obj, sum)

  message("Compute the number of clusters")
  track_nb_cluster_per_obj <-
    pbapply::pblapply(list_of_objects, function(object) {
      message(paste("Start object of class:", class(object), sep = " "))
      if (inherits(object, "phyloseq")) {
        ncol(object@otu_table)
      } else if (inherits(object, "matrix")) {
        ncol(object)
      } else if (inherits(object[[1]], "dada")) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (is.data.frame(object[[1]]) &&
        all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (inherits(object[[1]], "derep")) {
        length(unique(unlist(lapply(object, function(x) {
          names(x$uniques)
        }))))
      } else {
        NA
      }
    })

  message("Compute the number of samples")
  track_nb_sam_per_obj <-
    pbapply::pblapply(list_of_objects, function(object) {
      message(paste("Start object of class:", class(object), sep = " "))
      if (inherits(object, "phyloseq")) {
        nrow(object@otu_table)
      } else if (inherits(object, "matrix")) {
        nrow(object)
      } else if (inherits(object[[1]], "dada")) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (is.data.frame(object[[1]]) &&
        all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (inherits(object[[1]], "derep")) {
        length(object)
      } else if (is.character(object[1]) &&
        length(object[1]) == 1 &&
        file.exists(object[1])) {
        length(object)
      } else {
        NA
      }
    })


  track <- matrix(nrow = length(list_of_objects), ncol = 3)
  colnames(track) <- c(
    "nb_sequences",
    "nb_clusters",
    "nb_samples"
  )
  track <- as.data.frame(track)
  if (!is.null(obj_names)) {
    rownames(track) <- obj_names
  } else {
    rownames(track) <- names(list_of_objects)
  }

  track$nb_sequences <- unlist(track_nb_seq_per_obj)
  track$nb_clusters <- unlist(track_nb_cluster_per_obj)
  track$nb_samples <- unlist(track_nb_sam_per_obj)


  return(track)
}
################################################################################

################################################################################
#' Recluster sequences of an object of class `physeq`
#' (e.g. OTUs or ASV from dada)
#'
#' `r lifecycle::badge("maturing")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param method (default IdClusters)
#'   Set the clustering method.
#'   - `IdClusters` use the DECIPHER::IdClusters fonction,
#'   - `vsearch` use the vsearch software (https://github.com/torognes/vsearch/)
#'     with arguments `-cluster_fast` and `-strand both`
#' @param vsearchpath (default = "vsearch"): path to vsearch
#' @param id (default=0.97): level of identity to cluster
#' @param  tax_adjust: By default tax_adjust = 1L. See the man page
#'   of `speedyseq::merge_taxa_vec`. 
#' @details This function use the `speedyseq::merge_taxa_vec` function to
#'   merge taxa into clusters. By default tax_adjust = 1L. See the man page
#'   of `speedyseq::merge_taxa_vec`.
#'
#' @return A new object of class `physeq`
#'
#' @references
#'   VSEARCH can be downloaded from
#'   \url{https://github.com/torognes/vsearch}.
#'   More information in the associated publication
#'   \url{https://www.ncbi.nlm.nih.gov/pubmed/27781170}.
#'
#' @export

asv2otu <- function(physeq,
                    nproc = 1,
                    method = "IdClusters",
                    id = 0.97,
                    vsearchpath = "vsearch",
                    tax_adjust = 1
                    ) {
  dna <- Biostrings::DNAStringSet(physeq@refseq)

  if (method == "IdClusters") {
    ## Find clusters of ASVs to form the new OTUs
    clusters <- DECIPHER::IdClusters(dna,
      cutoff = 1 - id,
      # e.g. `cutoff = 0.03` for a 97% OTU
      processors = nproc
    )

    new_physeq <-
      speedyseq::merge_taxa_vec(physeq,
        clusters$cluster,
        tax_adjust = tax_adjust
      )
  } else if (method == "vsearch") {
    Biostrings::writeXStringSet(dna, "temp.fasta")

    cmd <-
      system2(
        vsearchpath,
        paste(
          " -cluster_fast temp.fasta -strand both",
          " -id ",
          id,
          " --centroids cluster.fasta",
          " --uc temp.uc",
          sep = ""
        ),
        stdout = T,
        stderr = T
      )
    vsearch_cluster_dna <-
      Biostrings::readDNAStringSet("cluster.fasta")

    pack_clusts <- utils::read.table("temp.uc", sep = "\t")
    colnames(pack_clusts) <-
      c(
        "type",
        "cluster",
        "width",
        "identity",
        "strand",
        "6",
        "7",
        "cigarAlignment",
        "query",
        "target"
      )

    clusters <- pack_clusts$cluster
    clusters <-
      tapply(
        clusters, paste(pack_clusts$cluster, pack_clusts$query),
        function(x) {
          x[1]
        }
      )

    new_physeq <-
      speedyseq::merge_taxa_vec(physeq,
        clusters,
        tax_adjust = tax_adjust
      )

    if (file.exists("temp.fasta")) {
      file.remove("temp.fasta")
    }
    if (file.exists("cluster.fasta")) {
      file.remove("cluster.fasta")
    }
    if (file.exists("temp.uc")) {
      file.remove("temp.uc")
    }
  }
  return(new_physeq)
}
################################################################################

################################################################################
#' Search for a list of sequence in a fasta file against physeq reference
#'   sequences using [vsearch](https://github.com/torognes/vsearch)
#'
#' `r lifecycle::badge("maturing")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param seq2search (required): path to fasta file
#' @param vsearchpath (default = "vsearch"): path to vsearch
#' @param id (default = 0.8): id for --usearch_global
#' @param iddef (default = 0): iddef for --usearch_global
#'
#' @return A dataframe with uc results (invisible)
#' @export

vsearch_search_global <- function(physeq,
                                  seq2search,
                                  vsearchpath = "vsearch",
                                  id = 0.8,
                                  iddef = 0) {
  dna <- Biostrings::DNAStringSet(physeq@refseq)

  Biostrings::writeXStringSet(dna, "temp.fasta")

  cmd <- system2(
    vsearchpath,
    paste(
      " --usearch_global ",
      here::here("temp.fasta"),
      " --db ",
      here::here(seq2search),
      " --uc",
      " temp.uc",
      " --id ",
      id,
      " --uc_allhits",
      " --strand both",
      " --iddef ",
      iddef,
      sep = ""
    )
  )

  pack_clusts <- utils::read.table("temp.uc", sep = "\t")
  colnames(pack_clusts) <- c(
    "type",
    "cluster",
    "width",
    "identity",
    "strand",
    "6",
    "7",
    "cigarAlignment",
    "query",
    "target"
  )

  return(invisible(pack_clusts))

  if (file.exists("temp.fasta")) {
    file.remove("temp.fasta")
  }
  if (file.exists("temp.uc")) {
    file.remove("temp.uc")
  }
}
################################################################################

################################################################################
#' Blast some sequence against `refseq` slot of a \code{\link{phyloseq-class}}
#'   object.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param seq2search (required): path to fasta file
#' @param blastpath (default = NULL): path to blast program
#' @param id_cut (default = 90): cut of in identity percent to keep result
#' @param bit_score_cut (default = 1e-10): cut of in bit score to keep result
#' @param unique_per_seq (default = FALSE) : only return the first match for
#'  each sequence in seq2search
#' @param score_filter (default = TRUE): does results are filter by score. If
#'   FALSE, `id_cut` and `bit_score_cut` are ignored
#' @param list_no_output_query (default to FALSE): does the result table include
#'   query sequences for which `blastn` does not find any correspondence
#'
#' @return  the blast table
#' @export
#'
#' @examples
#' \dontrun{
#' data(data_fungi)
#' blastpath <- "...YOUR_PATH_TO_BLAST..."
#' blast_to_phyloseq(data_fungi,
#'   seq2search = "inst/extdata/ex.fasta",
#'   blastpath = blastpath
#' )
#' }
blast_to_phyloseq <- function(physeq,
                              seq2search,
                              blastpath = NULL,
                              id_cut = 90,
                              bit_score_cut = 50,
                              unique_per_seq = FALSE,
                              score_filter = TRUE,
                              list_no_output_query = FALSE) {
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, "db.fasta")

  cmd <- system(paste(blastpath,
    "makeblastdb -dbtype nucl -in db.fasta -out dbase",
    sep = ""
  ))

  cmd2 <- system(
    paste(
      blastpath,
      "blastn -query ",
      seq2search,
      " -db dbase",
      " -out blast_result.txt",
      " -outfmt \"6 qseqid qlen sseqid slen",
      " length pident evalue bitscore qcovs\"",
      sep = ""
    )
  )

  if (file.info("blast_result.txt")$size > 0) {
    blast_tab <- utils::read.table(
      "blast_result.txt",
      sep = "\t",
      header = F,
      stringsAsFactors = F
    )
    file.remove("blast_result.txt")
    file.remove(list.files(pattern = "dbase"))
    file.remove("db.fasta")
  } else {
    file.remove("blast_result.txt")
    file.remove(list.files(pattern = "dbase"))
    file.remove("db.fasta")
    stop("None query sequences matched your phyloseq references sequences.")
  }


  names(blast_tab) <- c(
    "Query name",
    "Query seq. length",
    "Taxa name",
    "Taxa seq. length",
    "Alignment length",
    "% id. match",
    "e-value",
    "bit score",
    "Query cover"
  )

  blast_tab <- blast_tab[order(blast_tab[, "bit score"], decreasing = FALSE), ]

  if (unique_per_seq) {
    blast_tab <- blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if (score_filter) {
    blast_tab <- blast_tab[blast_tab[, "bit score"] > bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  if (list_no_output_query) {
    fastaFile <- Biostrings::readDNAStringSet(seq2search)
    seq_name <- names(fastaFile)
    no_output_query <- seq_name[!seq_name %in% blast_tab[1, ]]
    if (length(no_output_query) > 0) {
      mat_no_output_query <- matrix(NA,
        ncol = ncol(blast_tab),
        nrow = length(no_output_query)
      )
      mat_no_output_query[, 1] <- no_output_query
      colnames(mat_no_output_query) <- colnames(blast_tab)
      blast_tab <- rbind(blast_tab, mat_no_output_query)
    }
  } else {
    blast_tab <- blast_tab
  }

  return(blast_tab)
}
################################################################################

################################################################################
#' Save phyloseq object in the form of multiple csv tables.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param path (default: NULL) : a path to the folder to save the phyloseq object
#' @param Rdata (default: FALSE) : does the phyloseq object is also save in Rdata format
#'
#' @return One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#' and if present a phy_tree in Newick format
#' @export
#'
#' @examples
#' \dontrun{
#' write_phyloseq(data_fungi, path = "phyloseq")
#' }

write_phyloseq <- function(physeq, path = NULL, rdata = FALSE) {
  if (!dir.exists(path)){
    dir.create(file.path(path), recursive = TRUE)
  }
  if (!is.null(physeq@otu_table)) {
    utils::write.csv(physeq@otu_table, paste(path, "/otu_table.csv", sep = ""))
  }
  if (!is.null(physeq@refseq)) {
    utils::write.csv(physeq@refseq, paste(path, "/refseq.csv", sep = ""))
  }
  if (!is.null(physeq@tax_table)) {
    utils::write.csv(physeq@tax_table, paste(path, "/tax_table.csv", sep = ""))
  }
  if (!is.null(physeq@sam_data)) {
    utils::write.csv(physeq@sam_data, paste(path, "/sam_data.csv", sep = ""))
  }
  if (!is.null(physeq@phy_tree)) {
    ape::write.tree(physeq@phy_tree, paste(path, "/phy_tree.txt", sep = ""))
  }
  if (rdata){
    save(physeq, file=paste(path, "/physeq.RData", sep = ""))
  }
}
################################################################################

################################################################################
#' Read phyloseq object from multiple csv tables and a phylogenetic tree
#' in Newick format.
#'
#' `r lifecycle::badge("maturing")`
#'
#' @param path (required) : a path to the folder to read the phyloseq object
#' @param taxa_are_rows (required, default to FALSE) : see ?phyloseq for details
#'
#' @return One to four csv tables (refseq.csv, otu_table.csv, tax_table.csv, sam_data.csv)
#' and if present a phy_tree in Newick format. At least the otu_table.csv need to be present.
#' @export
#'
#' @examples
#' \dontrun{
#' read_phyloseq(path = "phyloseq_data")
#' }

read_phyloseq <- function(path = NULL, taxa_are_rows = FALSE) {

  if (file.exists(paste(path, "/otu_table.csv", sep = ""))) {
    otu_table_csv <- as.matrix(utils::read.csv(paste(path, "/otu_table.csv", sep = "")))
    rownames(otu_table_csv) <- otu_table_csv[,1]
    otu_table_csv <- otu_table_csv [,-1]
    otu_table_csv <- apply(otu_table_csv, 2, as.numeric)
    physeq <- phyloseq(otu_table(otu_table_csv, taxa_are_rows = taxa_are_rows))
  }
  if (file.exists(paste(path, "/refseq.csv", sep = ""))) {
    dna <- Biostrings::DNAStringSet(utils::read.csv2(paste(path, "/refseq.csv", sep = ""), sep=",")[,2])
    names(dna) <- utils::read.csv2(paste(path, "/refseq.csv", sep = ""), sep=",")[,1]
    physeq <- phyloseq::merge_phyloseq(physeq, refseq(dna))
  }
  if (file.exists(paste(path, "/tax_table.csv", sep = ""))) {
    tax_table_csv <- utils::read.csv(paste(path, "/tax_table.csv", sep = ""))
    rownames(tax_table_csv) <- tax_table_csv[,1]
    tax_table_csv <- as.matrix(tax_table_csv [,-1])
    physeq <- phyloseq::merge_phyloseq(physeq, tax_table(tax_table_csv) )
  }
   if (file.exists(paste(path, "/sam_data.csv", sep = ""))) {
    sam_data_csv <- utils::read.csv(paste(path, "/sam_data.csv", sep = ""))
    physeq <- phyloseq::merge_phyloseq(physeq, sample_data(sam_data_csv) )
  }

  if (!is.null(physeq@phy_tree)) {
    tree <- ape::read.tree(paste(path, "/phy_tree.txt", sep = ""))
    physeq <- phyloseq::merge_phyloseq(physeq, phy_tree(tree))
  }

  return(physeq)
}

################################################################################

################################################################################
#' Lulu reclustering of class `physeq`
#'
#' See https://www.nature.com/articles/s41467-017-01312-x for more information on the method.
#' `r lifecycle::badge("experimental")`
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for the clustering
#' @param id (default = 0.84): id for --usearch_global
#' @param vsearchpath (default = "vsearch"): path to vsearch
#' @param verbose (Logical): if true, print some additional messages
#' @param clean_physeq (Logical, default = FALSE) If true, empty samples and empty ASV are discarded before clustering.
#'
#' @return a list of for object
#' - "new_physeq": The new phyloseq object (class physeq)
#' - "discrepancy_vector": A vector of discrepancy showing for each taxonomic level the proportion of identic value
#'   before and after lulu reclustering. A value of 0.6 stands for 60% of ASV before re-clusering have
#'   identical value after re-clustering. In other word, 40% of ASV are assigned to a different taxonomic
#'   value. NA value are not counted as discrepancy.
#' - "res_lulu": A list of the result from the lulu function
#' - "merged_ASV": the data.frame used to merged ASV
#'
#' @export
#' @examples
#' \dontrun{
#' data(data_fungi_sp_known)
#' lulu_phyloseq(data_fungi_sp_known)
#' }
#' @author Adrien Taudière \email{adrien.taudiere@@zaclys.net}
#' @details
#' The version of LULU is a fork of Adrien Taudière (\url{https://github.com/adrientaudiere/lulu}) from \url{https://github.com/tobiasgf/lulu}
#' @references
#' - LULU : \url{https://github.com/adrientaudiere/lulu}
#'  forked from \url{https://github.com/tobiasgf/lulu}.
#' - VSEARCH can be downloaded from
#'  \url{https://github.com/torognes/vsearch}.

lulu_phyloseq <- function(physeq,
                          nproc = 1,
                          id = 0.84,
                          vsearchpath = "vsearch",
                          verbose = FALSE,
                          clean_physeq = FALSE) {
  if (clean_physeq) {
      physeq <- clean_physeq(physeq)
  }

  message("Start Vsearch usearch_global")
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, "temp.fasta")
  cmd <-
    system2(
      vsearchpath,
      paste(
        " --usearch_global temp.fasta --db temp.fasta --self --iddef 1",
        " -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10",
        " -id ",
        id,
        "  --userout match_list.txt",
        sep = ""
      ),
      stdout = T,
      stderr = T
    )

  match_list <- utils::read.csv(file = "match_list.txt", sep = "\t")

  if (!requireNamespace("lulu")) {
    requireNamespace(devtools)
    install_github("adrientaudiere/lulu")
  }

  message("Lulu algorithm")
  res_lulu <-
    lulu::lulu(data.frame(t(physeq@otu_table)), match_list)

  if (file.exists("temp.fasta")) {
    file.remove("temp.fasta")
  }
  if (file.exists("cluster.fasta")) {
    file.remove("cluster.fasta")
  }
  if (file.exists("temp.uc")) {
    file.remove("temp.uc")
  }
  if (file.exists("match_list.txt")) {
    file.remove("match_list.txt")
  }

  merged <- res_lulu$otu_map[res_lulu$otu_map$curated == "merged", ]
  merged <- merged[rownames(merged) != merged$parent_id, ]

  test_vector <- c()
  for (tax_rank in colnames(physeq@tax_table)) {
    test <- physeq@tax_table[rownames(merged), tax_rank] == physeq@tax_table[merged$parent_id, tax_rank]
    test_vector <- c(test_vector, sum(test, na.rm = T) / length(stats::na.exclude(test)))
  }

  names(test_vector) <- colnames(physeq@tax_table)

  new_physeq <- prune_taxa(taxa_names(physeq) %in% rownames(res_lulu$curated_table), physeq)
  new_physeq@otu_table <- otu_table(t(res_lulu$curated_table), taxa_are_rows = FALSE)
  sample_names(new_physeq) <- sample_names(physeq)

  if (verbose) {
    message(paste("The number of taxa decrease from ", ntaxa(physeq), " to ", ntaxa(new_physeq), ".", sep = ""))
    message("See the discrepancy_vector to verify the degree of discrepancy in taxonomy due to lulu re-clustering.")
  }
  return(list(
    "new_physeq" = new_physeq,
    "discrepancy_vector" = test_vector,
    "res_lulu" = res_lulu,
    "merged_ASV" = merged
  ))
}
################################################################################


