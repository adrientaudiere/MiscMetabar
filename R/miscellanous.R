################################################################################
#' Transform the otu_table of a \code{\link{phyloseq-class}} object into a
#'   \code{\link{phyloseq-class}} object with a binary otu_table.
#' @note  Useful to test if the results are not biased by sequences bias
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
#'
#' @inheritParams clean_pq
#' @param remove_space (logical; default TRUE): do we remove space?
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @author Adrien Taudière
#'
#' @return A  \code{\link{phyloseq-class}} object with simplified taxonomy
#' @export
simplify_taxo <- function(physeq, remove_space = TRUE) {
  taxo <- physeq@tax_table
  taxo <- gsub(".__", "", taxo, perl = TRUE)
  if (remove_space) {
    taxo <- gsub(" ", "", taxo)
    taxo <- gsub("\u00a0", "", taxo)
  }
  physeq@tax_table <- tax_table(taxo)
  return(physeq)
}
################################################################################

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


################################################################################
#' Convert a value (or a fraction x/y) in percentage
#'
#' @param x (required): value
#' @param y if y is set, compute the division of x by y
#' @param accuracy number of digits (number of digits after zero)
#' @param add_symbol if set to TRUE add the % symbol to the value
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
################################################################################


################################################################################
#' Count sequences in fasta or fastq file
#'
#' @description
#'  `r lifecycle::badge("experimental")`
#'   Use grep to count the number of line with only one '+' (fastq, fastq.gz)
#'   or lines starting with a '>' (fasta) to count sequences.
#'
#' @param file_path The path to a  fasta, fastq or fastq.gz file
#' @param folder_path The path to a folder with fasta, fastq or fastq.gz files
#' @param pattern A pattern to filter files in a folder. E.g. _R2_
#' @return the number of sequences
#' @author Adrien Taudière
#' @export
#' @examples
#' count_seq(file_path = system.file(
#'   "extdata",
#'   "ex.fasta",
#'   package = "MiscMetabar",
#'   mustWork = TRUE
#' ))
#' count_seq(folder_path = "extdata/", pattern = "*.fasta")
count_seq <- function(file_path = NULL, folder_path = NULL, pattern = NULL) {
  if (is.null(file_path) && is.null(folder_path)) {
    stop("You need to specify one of file_path or folder_path param!")
  } else if (!is.null(file_path) && !is.null(folder_path)) {
    stop("You need to specify either file_path or folder_path param not both!")
  } else if (!is.null(file_path) && is.null(folder_path)) {
    if (sum(get_file_extension(file_path) %in% "fasta") > 0) {
      seq_nb <- system(paste0("cat ", file_path, " | grep -ce '^>'"),
        intern = TRUE
      )
    } else if (sum(get_file_extension(file_path) %in% "fastq") > 0) {
      if (sum(get_file_extension(file_path) %in% "gz") > 0) {
        seq_nb <- system(paste0("zcat ", file_path, " | grep -ce '^+$'"),
          intern = TRUE
        )
      } else {
        seq_nb <- system(paste0("cat ", file_path, " | grep -ce '^+$'"),
          intern = TRUE
        )
      }
    } else {
      stop(paste0(
        "The file extension ",
        get_file_extension(file_path),
        " is not supported."
      ))
    }
  } else {
    seq_nb <- sapply(
      list.files(folder_path, full.names = TRUE, pattern = pattern),
      function(f) {
        count_seq(file_path = f)
      }
    )
  }
  return(as.numeric(seq_nb))
}

################################################################################


################################################################################
#' Funky palette color
#' @return a color palette
#' @param n a number of colors
#' @author Thibaut Jombart
#' @export
#'
funky_color <-
  grDevices::colorRampPalette(
    c(
      "#A6CEE3",
      "#1F78B4",
      "#B2DF8A",
      "#33A02C",
      "#FB9A99",
      "#E31A1C",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    )
  )
################################################################################


################################################################################

################################################################################
#' Subsample a fastq file copying the n_seq first sequences in a given folder
#'
#' @description
#'  `r lifecycle::badge("experimental")`
#'
#' @param fastq_files The path to one fastq file or a list of fastq files
#'   (see examples)
#' @param folder_output The path to a folder for output files
#' @param n_seq (int; default 1000) : Number of sequences kept (every sequence
#'   spread across 4 lines)
#' @return Nothing, create subsampled fastq files in a folder
#' @author Adrien Taudière
#' @export
#' @examples
#' \dontrun{
#' ex_file <- system.file("extdata", "ex_R1_001.fastq.gz",
#'   package = "MiscMetabar",
#'   mustWork = TRUE
#' )
#' subsample_fastq(ex_file, "your_path_to_output")
#' subsample_fastq(list_fastq_files("extdata"), "your_path_to_output", n = 10)
#' }
subsample_fastq <- function(fastq_files,
                            folder_output = "subsample",
                            n_seq = 1000) {
  for (f in unlist(fastq_files)) {
    if (!dir.exists(folder_output)) {
      dir.create(folder_output)
    }
    writeLines(readLines(f, n = n_seq * 4), con = paste0(
      folder_output, "/",
      basename(f)
    ))
  }
}

################################################################################

################################################################################
#' Install a package if not present
#'
#' @description
#'  `r lifecycle::badge("experimental")`
#'
#' @param pkg The name of the package
#' @param use_pak (logical, default TRUE) Use of `pak::pkg_install()`. If FALSE
#'   use the base `install.package()` function or the function
#'   `BiocManager::install()` if bioconductor_pkg is true or the function
#' @param bioconductor_pkg (logical, default FALSE). If use_pak is TRUE,
#'   do nothing, else use `BiocManager::install()` to install the package.
#' @param github_pkg (logical, default FALSE). If use_pak is TRUE,
#'   do nothing, else use `devtools::install_github` to install the package.
#' @param verbose (logical, default FALSE) Does the function print message?
#' @return Nothing
#' @author Adrien Taudière
#' @export
#' @examples
#' \dontrun{
#' install_pkg_needed("ggVennDiagram")
#' }
install_pkg_needed <- function(pkg, use_pak = TRUE, bioconductor_pkg = FALSE,
                               github_pkg = FALSE, verbose = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (verbose) {
      message(paste0("Installation of the package : ", pkg))
    }
    if (use_pak) {
      if (!requireNamespace("pak", quietly = TRUE)) {
        install.packages("pak")
      }
      pak::pkg_install(pkg)
    } else {
      if (bioconductor_pkg) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else if (github_pkg) {
        if (!requireNamespace("devtools", quietly = TRUE)) {
          install.packages("devtools", quietly = TRUE)
        }
        devtools::install_github(pkg)
      } else {
        install.packages(pkg)
      }
    }
  } else {
    if (verbose) {
      message(paste0(pkg, " is already present."))
    }
  }
}
