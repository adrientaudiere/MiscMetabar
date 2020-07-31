################################################################################
#' Add dna in `refseq` slot of a physeq object using taxa names and renames taxa
#'   using ASV_1, ASV_2, …
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object obtained
#'   using the `dada2` package
#'
#' @return A new \code{\link{phyloseq-class}} object with `refseq` slot and new
#'   taxa names
#' @export

add_dna_to_phyloseq <- function(physeq) {
  if (!validObject(physeq) | class(physeq) != "phyloseq") {
    stop("The physeq argument is not a valid phyloseq object.")
  }
  else {
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
#' Track the number of reads (= sequences), samples and cluster (e.g. ASV)
#' from various objects including dada-class and derep-class.
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
#'
#' @return The number of sequences, clusters (e.g. OTUs, ASVs) and samples for
#'   each object.
#' @export

track_wkflow <- function(list_of_objects, obj_names = NULL) {
  message("Compute the number of sequences")
  if (!is.null(obj_names)) {
    names(list_of_objects) <- obj_names
  }
  track_nb_seq_per_obj <-
    pbapply::pblapply(list_of_objects, function(object) {
      message(paste("Start object of class:", class(object), sep = " "))
      if (class(object) == "phyloseq") {
        rowSums(object@otu_table)
      } else if (class(object) == "matrix") {
        sum(object, na.rm = TRUE)
      } else if (is.character(object[1]) &&
                 length(object[1]) == 1 &&
                 file.exists(object[1])) {
        if (summary(file(object[[1]]))$class == "gzfile") {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("zcat ", x, " | grep -c '^+$'", sep = ""),
                              intern = TRUE))
          })
        } else if (grepl("\\.fastq$", object[1])) {
          pbapply::pbsapply(object, function(x) {
            as.numeric(system(paste("cat ", x, " | grep -c '^+$'", sep = ""),
                              intern = TRUE))
          })
        } else {
          stop("Files must be either gzfile or .fastq")
        }
      }
      else {
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
      if (class(object) == "phyloseq") {
        ncol(object@otu_table)
      } else if (class(object) == "matrix") {
        ncol(object)
      } else if (class(object[[1]]) == "dada") {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (is.data.frame(object[[1]]) &&
                 all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[2]
      } else if (class(object[[1]]) == "derep") {
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
      if (class(object) == "phyloseq") {
        nrow(object@otu_table)
      } else if (class(object) == "matrix") {
        nrow(object)
      } else if (class(object[[1]]) == "dada") {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (is.data.frame(object[[1]]) &&
                 all(c("sequence", "abundance") %in% colnames(object[[1]]))) {
        dim(suppressMessages(dada2::makeSequenceTable(object)))[1]
      } else if (class(object[[1]]) == "derep") {
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
  colnames(track) <- c("nb_sequences",
                       "nb_clusters",
                       "nb_samples")
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
                    vsearchpath = "vsearch") {
  dna <- Biostrings::DNAStringSet(physeq@refseq)

  if (method == "IdClusters") {
    ## Find clusters of ASVs to form the new OTUs
    aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
    d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
    clusters <- DECIPHER::IdClusters(d,
                                     method = "complete",
                                     cutoff = 1 - id,
                                     # e.g. `cutoff = 0.03` for a 97% OTU
                                     processors = nproc)

    new_phyloseq <-
      speedyseq::merge_taxa_vec(physeq,
                                clusters$cluster,
                                tax_adjust = 1)
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
      tapply(clusters, paste(pack_clusts$cluster, pack_clusts$query),
             function(x) {
               x[1]
             })

    new_phyloseq <-
      speedyseq::merge_taxa_vec(physeq,
                                clusters,
                                tax_adjust = 1)
    new_phyloseq@refseq <- vsearch_cluster_dna

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
  return(new_phyloseq)
}
################################################################################

################################################################################
#' Search for a list of sequence in a fasta file against physeq reference
#'   sequences using [vsearch](https://github.com/torognes/vsearch)
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
#' @return
#' @export
#'
#' @examples
#' data(data_fungi)
#' blastpath =
#' "/home/adrien/Bureau/ncbi-blast-2.10.0+-x64-linux/ncbi-blast-2.10.0+/bin/"
#' blast_to_phyloseq(data_fungi, seq2search = "data/ex.fasta",
#'                  blastpath = blastpath)
blast_to_phyloseq <- function(physeq,
                              seq2search,
                              blastpath = NULL,
                              id_cut = 90,
                              bit_score_cut = 1e-10,
                              unique_per_seq = FALSE,
                              score_filter = TRUE,
                              list_no_output_query = FALSE) {
  dna <- Biostrings::DNAStringSet(physeq@refseq)
  Biostrings::writeXStringSet(dna, "db.fasta")

  cmd <- system(paste(blastpath,
                       "makeblastdb -dbtype nucl -in db.fasta -out dbase",
                       sep = "")
                 )

  cmd2 <- system(
    paste(
      blastpath,
      "blastn -query ",
      seq2search,
      " -db dbase",
      " -out blast_result.txt",
      " -outfmt \"6 qseqid qlen sseqid slen length pident evalue bitscore qcovs\"",
      sep = "")
  )

  if (file.info("blast_result.txt")$size > 0){
    blast_tab <- read.table(
      "blast_result.txt",
      sep = "\t",
      header = F,
      stringsAsFactors = F
    )
  } else {
    stop("None query sequences matched your phyloseq references sequences.")
  }

  file.remove("blast_result.txt")
  file.remove(list.files(pattern = "dbase"))
  file.remove("db.fasta")

  names(blast_tab) <- c("Query name",
                        "Query seq. length",
                        "Taxa name",
                        "Taxa seq. length",
                        "Alignment length",
                        "% id. match",
                        "bit score",
                        "e-value",
                        "Query cover"
                        )

  blast_tab <- blast_tab[order(blast_tab[,"bit score"], decreasing = FALSE), ]

  if (unique_per_seq) {
    blast_tab <-  blast_tab[which(!duplicated(blast_tab[, 1])), ]
  }

  if(score_filter){
    blast_tab <- blast_tab[blast_tab[, "bit score"] < bit_score_cut, ]
    blast_tab <- blast_tab[blast_tab[, "% id. match"] > id_cut, ]
  } else {
    blast_tab <- blast_tab
  }

  if (list_no_output_query) {
    fastaFile <- Biostrings::readDNAStringSet(seq2search)
    seq_name <- names(fastaFile)
    no_output_query <- seq_name[!seq_name %in% blast_tab[1, ]]
    if(length(no_output_query)>0) {
      mat_no_output_query <- matrix(NA,
                                    ncol = ncol(blast_tab),
                                    nrow = length(no_output_query))
      mat_no_output_query[,1] <- no_output_query
      colnames(mat_no_output_query) <- colnames(blast_tab)
      blast_tab <- rbind(blast_tab, mat_no_output_query)
    }
  } else {
    blast_tab <- blast_tab
  }
}
################################################################################

################################################################################
#' Save phyloseq object in the form of multiple csv tables. EXPERIMENTAL
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param path (defaut: NULL) : a path to the folder to save the phyloseq object
#'
#' @return
#' @export
#'
#' @examples
write_phyloseq <- function(physeq, path = NULL) {
  if (!is.null(physeq@otu_table)) {
    write.csv(physeq@otu_table, paste(path, "otu_table.csv", sep=""))
  }
  if (!is.null(physeq@refseq)) {
    write.csv(physeq@refseq, paste(path, "refseq.csv", sep=""))
  }
  if (!is.null(physeq@tax_table)) {
    write.csv(physeq@tax_table, paste(path, "tax_table.csv", sep=""))
  }
  if (!is.null(physeq@sam_data)) {
    write.csv(physeq@sam_data, paste(path, "sam_data.csv", sep=""))
  }
}
################################################################################


