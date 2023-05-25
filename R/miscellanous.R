################################################################################
#' Transform the otu_table of a \code{\link{phyloseq-class}} object into a
#'   \code{\link{phyloseq-class}} object with à binary otu_table.
#' @note  Useful to test if the results are not biaised by sequences bias
#'   that appended during PCR or NGS pipeline.
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @inheritParams clean_pq
#' @param min_number (int) the minimum number of sequences to put
#'   a 1 in the otu table.
#' @author Adrien Taudière
#'
#' @return A \code{physeq} object with only 0/1 in the OTU table
#' @export
#' @examples
#' data(enterotype)
#' enterotype_bin <- as_binary_otu_table(enterotype)
as_binary_otu_table <- function(physeq, min_number = 1) {
  if (!inherits(physeq, "phyloseq")) {
    stop("physeq must be a phyloseq object")
  }
  res <- physeq
  res@otu_table[res@otu_table >= min_number] <- 1
  res@otu_table[res@otu_table < min_number] <- 0
  return(res)
}
################################################################################


################################################################################
#' Compute paired distances among matrix (e.g. otu_table)
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @note the first column of the first matrix is compare to the first column of
#'   the second matrix, the second column of the first matrix is compare to the
#'   second column of the second matrix and so on.
#' @param x (required) A first matrix.
#' @param y (required) A second matrix.
#' @param method (default: 'bray') the method to use internally in the vegdist
#'   function.
#' @param nperm (default: 99) The number of permutations
#' @param ... others argument for `vegan::vegdist` function
#'
#' @author Adrien Taudière
#'
#' @return A list of length two : (i) a vector of observed distance ($obs) and
#'   (ii) a matrix of the distance after randomization ($null)
#' @export
#' @seealso \code{\link[vegan]{vegdist}}

dist_bycol <- function(x,
                       y,
                       method = "bray",
                       nperm = 99,
                       ...) {
  x <- as.matrix(unclass(x))
  y <- as.matrix(unclass(y))

  if (nrow(x) != nrow(y) ||
    ncol(x) != ncol(y)) {
    stop("x and y must be of the same dimension")
  }

  res <- list()
  res$obs <- rep(NA, ncol(x))
  res$null <- list(length = nperm)

  for (i in seq_len(ncol(x))) {
    res$obs[i] <-
      vegan::vegdist(rbind(x[, i], y[, i]), method = method, ...)
  }

  for (n in 1:nperm) {
    y_null <- y[, sample(seq_len(ncol(y)), replace = FALSE)]
    res$null[[n]] <- rep(NA, ncol(x))
    for (i in seq_len(ncol(x))) {
      res$null[[n]][i] <-
        vegan::vegdist(rbind(x[, i], y_null[, i]), method = method, ...)
    }
    message(n)
  }

  names(res$obs) <- colnames(x)
  return(res)
}
################################################################################


################################################################################
#' List the size of all objects of the GlobalEnv.
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Code from https://tolstoy.newcastle.edu.au/R/e6/help/09/01/1121.html
#'
#' @aliases all_object_size
#' @return a list of size
#' @export
all_object_size <- function() {
  return(sort(sapply(ls(envir = .GlobalEnv), function(x) {
    utils::object.size(get(x))
  })))
}
################################################################################



################################################################################
#' Simplify taxonomy by removing some unused characters such as "k__"
#' @inheritParams clean_pq
#'
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @author Adrien Taudière
#'
#' @return A  \code{\link{phyloseq-class}} object with simplified taxonomy
#' @export
simplify_taxo <- function(physeq) {
  taxo <- physeq@tax_table
  taxo <- gsub(".__", "", taxo, perl = TRUE)
  physeq@tax_table <- taxo
  return(physeq)
}

################################################################################
#' Get the extension of a file
#'
#' @param file (required): path to a file
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @author Adrien Taudière
#'
#' @return A  \code{\link{phyloseq-class}} object with simplified taxonomy
#' @export
get_file_extension <- function(file) {
  file_ext <- strsplit(basename(file), ".", fixed = TRUE)[[1]][-1]
  return(file_ext)
}

################################################################################
#' Convert a value (or a fraction x/y) in percentage
#'
#' @param x (required): value
#' @param y if y is set, compute the division of x by y
#' @param accuracy: number of digits (number of digits after zero)
#' @param add_symbol: if set to TRUE add the % symbol to the value
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @author Adrien Taudière
#'
#' @return The percentage value (number or character if add_symbol
#'   is set to TRUE)
#' @export
perc <- function(x, y = NULL, accuracy = 0, add_symbol = FALSE) {
  if (is.null(y)) {
    res <- round(x * 100, digits = accuracy)
  } else {
    res <- round(x / y * 100, digits = accuracy)
  }

  if (add_symbol) {
    res <- paste0(res, "%")
  }
  return(res)
}

#' Count sequences in fasta or fastq file
#'
#' @description
#'  `r lifecycle::badge("experimental")`
#'   Use grep to count the number of line with only one '+' (fastq, fastq.gz)
#'   or lines starting with a '>' (fasta) to count sequences.
#'
#' @param file_path The path to a  fasta, fastq or fastq.gz file
#' @param folder_path The path to a folder with fasta, fastq or fastq.gz files
#' 
#' @return the number of sequences
#' @author Adrien Taudière
#' @export
#' @examples 
#' count_seq(file = "inst/extdata/ex.fasta")
#' count_seq(folder = "inst/extdata/")
count_seq <- function(file_path = NULL, folder_path = NULL) {
  if(is.null(file_path) && is.null(folder_path)){
    stop("You need to specify one of file_path or folder_path param!")
  } else if(!is.null(file_path) && !is.null(folder_path)){
    stop("You need to specify either file_path or folder_path param not both!")
  } else if(!is.null(file_path) && is.null(folder_path)){
    if (sum(get_file_extension(file_path) %in% "fasta")>0) {
      seq_nb <- system(paste0("cat ", file_path, " | grep -ce '^>'"),
        intern = TRUE
      )
    } else if (sum(get_file_extension(file_path) %in% "fastq")>0) {
      if (sum(get_file_extension(file_path) %in% "gz")>0) {
        seq_nb <- system(paste0("zcat ", file_path, " | grep -ce '^+$'"),
          intern = TRUE
        )
      } else {
        seq_nb <- system(paste0("cat ", file_path, " | grep -ce '^+$'"),
          intern = TRUE
        )
      }
    } else {
      stop(paste0("The file extension", 
                  get_file_extension(file_path), 
                  "is not supported."))
    }
  } else {
    seq_nb <- sapply(list.files(folder_path, full.names = TRUE), 
      function(f) {
        count_seq(file_path = f)
      }
    )
  }
  return(seq_nb)
}
