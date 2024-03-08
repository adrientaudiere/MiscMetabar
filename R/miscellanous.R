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
#' @importFrom utils object.size
#' @aliases all_object_size
#' @return a list of size
#' @export
all_object_size <- function() {
  return(sort(vapply(ls(envir = .GlobalEnv), function(x) {
    utils::object.size(get(x))
  }, numeric(1))))
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
#' @param file_path (required): path to a file
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' @author Adrien Taudière
#'
#' @return A  \code{\link{phyloseq-class}} object with simplified taxonomy
#' @export
get_file_extension <- function(file_path) {
  file_ext <- strsplit(basename(file_path), ".", fixed = TRUE)[[1]][-1]
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
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' count_seq(file_path = system.file(
#'   "extdata",
#'   "ex.fasta",
#'   package = "MiscMetabar",
#'   mustWork = TRUE
#' ))
#' count_seq(
#'   folder_path = system.file("extdata", package = "MiscMetabar"),
#'   pattern = "*.fasta"
#' )
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
    seq_nb <- vapply(
      list.files(folder_path, full.names = TRUE, pattern = pattern),
      function(f) {
        count_seq(file_path = f)
      }, numeric(1)
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
#' Subsample a fastq file copying the n_seq first sequences in a given folder
#'
#' @description
#'  `r lifecycle::badge("experimental")`
#'
#' @param fastq_files The path to one fastq file or a list of fastq files
#'   (see examples)
#' @param folder_output The path to a folder for output files
#' @param nb_seq (int; default 1000) : Number of sequences kept (every sequence
#'   spread across 4 lines)
#' @return Nothing, create subsampled fastq files in a folder
#' @author Adrien Taudière
#' @export
#' @examples
#' \donttest{
#' ex_file <- system.file("extdata", "ex_R1_001.fastq.gz",
#'   package = "MiscMetabar",
#'   mustWork = TRUE
#' )
#' subsample_fastq(ex_file, paste0(tempdir(), "/output_fastq"))
#' subsample_fastq(list_fastq_files("extdata"), paste0(tempdir(), "/output_fastq"), n = 10)
#' unlink(paste0(tempdir(), "/output_fastq"), recursive = TRUE)
#' }
subsample_fastq <- function(fastq_files,
                            folder_output = "subsample",
                            nb_seq = 1000) {
  for (f in unlist(fastq_files)) {
    if (!dir.exists(folder_output)) {
      dir.create(folder_output)
    }
    writeLines(readLines(f, n = nb_seq * 4), con = paste0(
      folder_output, "/",
      basename(f)
    ))
  }
}

################################################################################

################################################################################
#' Test if cutadapt is installed.
#'
#' @description
#'  `r lifecycle::badge("maturing")`
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param args_before_cutadapt : (String) A one line bash command to run before
#' to run cutadapt. For examples, "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv &&" allow to bypass the conda init which asks to restart the shell
#' @export
#' @return A logical that say if cutadapt is install in
#'
#' @examples
#' MiscMetabar::is_cutadapt_installed()
#' @author Adrien Taudière

is_cutadapt_installed <- function(args_before_cutadapt = "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv && ") {
  writeLines(paste0(args_before_cutadapt, " cutadapt -h"), paste0(tempdir(), "/script_cutadapt.sh"))
  cutadapt_error_or_not <- try(system(paste0("bash ", tempdir(), "/script_cutadapt.sh 2>&1"), intern = TRUE), silent = T)
  unlink(paste0(tempdir(), "/script_cutadapt.sh"))

  return(class(cutadapt_error_or_not) != "try-error")
}

#' Test if falco is installed.
#'
#' @description
#'  `r lifecycle::badge("maturing")`
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: falco) Path to falco
#' @export
#' @return A logical that say if falco is install in
#'
#' @examples
#' MiscMetabar::is_falco_installed()
#' @author Adrien Taudière

is_falco_installed <- function(path = "falco") {
  return(class(try(system(paste0(path, " 2>&1"), intern = TRUE),
    silent = TRUE
  )) != "try-error")
}

#' Test if swarm is installed.
#'
#' @description
#'  `r lifecycle::badge("maturing")`
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: swarm) Path to falco
#' @export
#' @return A logical that say if swarm is install in
#'
#' @examples
#' MiscMetabar::is_swarm_installed()
#' @author Adrien Taudière

is_swarm_installed <- function(path = "swarm") {
  return(class(try(system(paste0(path, " -h 2>&1"), intern = TRUE),
    silent = TRUE
  )) != "try-error")
}

#' Test if vsearch is installed.
#'
#' @description
#'  `r lifecycle::badge("maturing")`
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: vsearch) Path to vsearch
#' @export
#' @return A logical that say if vsearch is install in
#'
#' @examples
#' MiscMetabar::is_vsearch_installed()
#' @author Adrien Taudière

is_vsearch_installed <- function(path = "vsearch") {
  return(class(try(system(paste0(path, " 2>&1"), intern = TRUE),
    silent = TRUE
  )) != "try-error")
}

#' Test if mumu is installed.
#'
#' @description
#'  `r lifecycle::badge("maturing")`
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: mumu) Path to mumu
#' @export
#' @return A logical that say if mumu is install in
#'
#' @examples
#' MiscMetabar::is_mumu_installed()
#' @author Adrien Taudière

is_mumu_installed <- function(path = "mumu") {
  return(class(try(system(paste0(path, " 2>&1"), intern = TRUE),
    silent = TRUE
  )) != "try-error")
}


#' Test if krona is installed.
#'
#' @description
#'  `r lifecycle::badge("maturing")`
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: krona) Path to krona
#' @export
#' @return A logical that say if krona is install in
#'
#' @examples
#' MiscMetabar::is_krona_installed()
#' @author Adrien Taudière

is_krona_installed <- function(path = "ktImportKrona") {
  return(class(try(system(paste0(path, " 2>&1"), intern = TRUE),
    silent = TRUE
  )) != "try-error")
}
################################################################################
