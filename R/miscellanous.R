################################################################################
#' Transform the otu_table of a \code{\link{phyloseq-class}} object into a
#'   \code{\link{phyloseq-class}} object with à binary otu_table.
#' @note  Usefull to test if the results are not biaised by sequences bias
#'   that appended during PCR or NGS pipeline.
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param minNumber (default = 1): the minimum number of sequences to put
#' a 1 in the otu table.
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
#' @note the first column of the first matrix is compare to the first column of
#'   the second matrix, the second column of the first matrix is compare to the
#'   second column of the second matrix and so on.
#' @param x (required): A first matrix.
#' @param y (required): A second matrix.
#' @param method (default = 'bray'): the method to use internally in the vegdist
#'   function.
#' @param nperm (default = 99): The number of permutations
#' @param ... (optional): others argument for `vegan::vegdist` function
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

  if (nrow(x) != nrow(y) |
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
    print(n)
  }

  names(res$obs) <- colnames(x)
  return(res)
}
################################################################################


################################################################################
#' List the size of all object of the GlobalEnv.
#' Code from https://tolstoy.newcastle.edu.au/R/e6/help/09/01/1121.html
#' @aliases all_object_size
#' @return a list of size
#' @export
all_object_size <- function() {
  return(sort(sapply(ls(envir = .GlobalEnv), function(x) {
    object.size(get(x))
  })))
}
################################################################################
