################################################################################
#' Transform the otu_table of a \code{\link[phyloseq]{phyloseq-class}} object into a
#'   \code{\link[phyloseq]{phyloseq-class}} object with a binary otu_table.
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' Useful to test if the results are not biased by sequences bias
#'   that appended during PCR or NGS pipeline.
#'
#' @inheritParams clean_pq
#' @param min_number (int) the minimum number of sequences to put
#'   a 1 in the OTU table.
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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' May be used to verify ecological distance among samples.
#'
#' @note the first column of the first matrix is compare to the first column of
#'   the second matrix, the second column of the first matrix is compare to the
#'   second column of the second matrix and so on.
#' @param x (required) A first matrix.
#' @param y (required) A second matrix.
#' @param method (default: 'bray') the method to use internally in the vegdist
#'   function.
#' @param nperm (int) The number of permutations to perform.
#' @param ... Additional arguments passed on to`vegan::vegdist` function
#'
#' @author Adrien Taudière
#'
#' @return A list of length two : (i) a vector of observed distance ($obs) and
#'   (ii) a matrix of the distance after randomization ($null)
#' @export
#' @seealso \code{\link[vegan]{vegdist}}
#' @examplesIf rlang::is_installed("vegan")
#' m1 <- matrix(runif(9), nrow = 3)
#' m2 <- matrix(runif(9), nrow = 3)
#' dist_bycol(m1, m2, nperm = 9)
dist_bycol <- function(x, y, method = "bray", nperm = 99, ...) {
  x <- as.matrix(unclass(x))
  y <- as.matrix(unclass(y))

  if (
    nrow(x) != nrow(y) ||
      ncol(x) != ncol(y)
  ) {
    stop("x and y must be of the same dimension")
  }

  res <- list(obs = rep(NA, ncol(x)), null = list(length = nperm))

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
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' Code from https://tolstoy.newcastle.edu.au/R/e6/help/09/01/1121.html
#'
#' @importFrom utils object.size
#' @aliases all_object_size
#' @return a list of size
#' @export
#' @examples
#' all_object_size()
all_object_size <- function() {
  return(sort(vapply(
    ls(envir = .GlobalEnv),
    function(x) {
      utils::object.size(get(x))
    },
    numeric(1)
  )))
}
################################################################################

################################################################################
#' Simplify taxonomy by removing some unused characters such as "k__"
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Internally used in [clean_pq()]
#' @inheritParams clean_pq
#' @param pattern_to_remove (a vector of character) regex patterns passed to
#'   [base::gsub()]: the matched *substring* is deleted from the cell value and
#'   the rest of the string is kept (e.g. `".__"` turns `"k__Fungi"` into
#'   `"Fungi"`).
#' @param ranks_for_pattern_to_remove (character vector or NULL; default all
#'   ranks) column names in `tax_table` to which `pattern_to_remove` is
#'   applied. Pass `NULL` to skip this operation on all columns.
#' @param ranks_to_remove_space (character vector or NULL; default all ranks
#'   whose name does not contain `"Species"`) column names from which ASCII
#'   spaces and non-breaking spaces (U+00A0) are stripped. Pass `NULL` to
#'   skip space removal entirely.
#' @param ranks_to_remove_NA (character vector or NULL; default all ranks)
#'   column names from which the literal string `"NA"` (case-sensitive) is
#'   removed. Pass `NULL` to skip this operation. **Breaking change from
#'   v0.16:** the old `remove_NA = FALSE` default is now `ranks_to_remove_NA`
#'   defaulting to all ranks; pass `NULL` to reproduce the old behaviour.
#' @param pattern_to_NA (character; default NULL): a regex; if an entire cell
#'   value matches, the *whole cell* is replaced with `NA` (nothing from the
#'   original value is kept). Designed for PR2-style placeholder unknowns such
#'   as `Embryophyceae_X`, `Embryophyceae_XX`, `Embryophyceae_XXX`,
#'   `Embryophyceae_XXX_sp.`, or `Mortierella_sp.`. Use `"_X+$|_sp\\.$"` to
#'   cover all such patterns: `_X+$` catches rank-filler X's; `_sp\\.$` catches
#'   any genus-only species placeholder.
#' @param ranks_for_pattern_to_NA (character vector or NULL; default all ranks)
#'   column names to which `pattern_to_NA` is applied. Pass `NULL` to skip
#'   this operation on all columns.
#' @author Adrien Taudière
#'
#' @return A  \code{\link[phyloseq]{phyloseq-class}} object with simplified taxonomy
#' @export
#' @examples
#' d_fm <- data_fungi_mini
#' d_fm@tax_table[, "Species"] <- paste0(rep(
#'   c("s__", "s:"),
#'   ntaxa(d_fm) / 2
#' ), d_fm@tax_table[, "Species"])
#'
#' # First column is the new vector of Species,
#' # second column is the column before simplification
#' cbind(
#'   simplify_taxo(d_fm)@tax_table[, "Species"],
#'   d_fm@tax_table[, "Species"]
#' )
#' # Apply pattern_to_remove only to Genus and Species columns
#' cbind(
#'   simplify_taxo(d_fm,
#'     ranks_for_pattern_to_remove = c("Genus", "Species")
#'   )@tax_table[, "Species"],
#'   d_fm@tax_table[, "Species"]
#' )
#' \dontrun{
#' # Replace PR2 placeholder unknowns (_X, _XX, _XXX, _XXX_sp., Genus_sp.) with NA
#' simplify_taxo(pq_pr2, pattern_to_NA = "_X+$|_sp\\.$")
#' # Apply pattern_to_NA only to the Species column
#' simplify_taxo(pq_pr2,
#'   pattern_to_NA = "_X+$|_sp\\.$",
#'   ranks_for_pattern_to_NA = "Species"
#' )
#' }
simplify_taxo <- function(
  physeq,
  pattern_to_remove = c(".__", ".*:"),
  ranks_for_pattern_to_remove = phyloseq::rank_names(physeq),
  ranks_to_remove_space = phyloseq::rank_names(physeq)[
    !grepl("Species", phyloseq::rank_names(physeq))
  ],
  ranks_to_remove_NA = phyloseq::rank_names(physeq),
  pattern_to_NA = NULL,
  ranks_for_pattern_to_NA = phyloseq::rank_names(physeq)
) {
  taxo <- as(physeq@tax_table, "matrix")

  if (!is.null(ranks_for_pattern_to_remove)) {
    for (p in pattern_to_remove) {
      for (col in ranks_for_pattern_to_remove) {
        taxo[, col] <- gsub(p, "", taxo[, col])
      }
    }
  }

  if (!is.null(ranks_to_remove_space)) {
    for (col in ranks_to_remove_space) {
      taxo[, col] <- gsub(" ", "", taxo[, col])
      taxo[, col] <- gsub("\u00a0", "", taxo[, col])
    }
  }

  if (!is.null(ranks_to_remove_NA)) {
    for (col in ranks_to_remove_NA) {
      taxo[, col] <- gsub("NA", "", taxo[, col], ignore.case = FALSE)
    }
  }

  if (!is.null(pattern_to_NA) && !is.null(ranks_for_pattern_to_NA)) {
    for (col in ranks_for_pattern_to_NA) {
      taxo[grepl(pattern_to_NA, taxo[, col]), col] <- NA
    }
  }

  physeq@tax_table <- tax_table(taxo)
  return(physeq)
}
################################################################################

################################################################################
#' Get the extension of a file
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' Internally used in [count_seq()].
#' Warning: don't work when there is '.' in the name of the
#'   file before the extension
#' @param file_path (required) path to a file
#' @author Adrien Taudière
#'
#' @return The extension of a file.
#' @export
#' @examples
#' get_file_extension("myfile.fasta")
get_file_extension <- function(file_path) {
  if (stringr::str_count(file_path, "\\.") == 0) {
    stop("There is no '.' inside your file path: ", file_path)
  }
  if (stringr::str_count(file_path, "\\.") > 1) {
    warning("There is more than one '.' inside your file path: ", file_path)
  }
  file_ext <- strsplit(basename(file_path), ".", fixed = TRUE)[[1]][-1]
  return(file_ext)
}
################################################################################

################################################################################
#' Convert a value (or a fraction x/y) in percentage
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Mostly for internal use.
#'
#' @param x (required) value
#' @param y if y is set, compute the division of x by y
#' @param accuracy number of digits (number of digits after zero)
#' @param add_symbol if set to TRUE add the % symbol to the value
#'
#' @author Adrien Taudière
#'
#' @return The percentage value (number or character if add_symbol
#'   is set to TRUE)
#' @export
#' @examples
#' perc(0.75)
#' perc(3, 10)
#' perc(0.75, add_symbol = TRUE)
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
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Use grep to count the number of line with only one '+' (fastq, fastq.gz)
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
    if (sum(get_file_extension(file_path) == "fasta") > 0) {
      if (sum(get_file_extension(file_path) == "gz") > 0) {
        seq_nb <- system(
          paste0("zcat ", file_path, " | grep -ce '^>'"),
          intern = TRUE
        )
      } else {
        seq_nb <- system(
          paste0("cat ", file_path, " | grep -ce '^>'"),
          intern = TRUE
        )
      }
    } else if (sum(get_file_extension(file_path) == "fastq") > 0) {
      if (sum(get_file_extension(file_path) == "gz") > 0) {
        seq_nb <- system(
          paste0("zcat ", file_path, " | grep -ce '^+$'"),
          intern = TRUE
        )
      } else {
        seq_nb <- system(
          paste0("cat ", file_path, " | grep -ce '^+$'"),
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
      },
      numeric(1)
    )
  }
  return(as.numeric(seq_nb))
}

################################################################################

################################################################################
#' Funky palette color
#' @return a color palette
#' @param n a number of colors
#' @author Thibaut Jombart in `adegenet` package
#' @export
#' @seealso The R package RColorBrewer, proposing a nice selection of color palettes. The viridis package, with many excellent palettes
#' @examples
#' funky_color(5)
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
#' Translates a factor into colors.
#' @param x	a numeric vector (for num2col) or a vector converted to a factor (for fac2col).
#' @param col.pal (default funky_color)	a function generating colors according to a given palette.
#' @param na.col (default grey) the color to be used for missing values (NAs)
#' @param seed (default NULL) a seed for R's random number generated, used to fix the random permutation of colors in the palette used; if NULL, no randomization is used and the colors are taken from the palette according to the ordering of the levels
#' @return a color vector
#' @author Thibaut Jombart in `adegenet` package
#' @export
#' @seealso The R package RColorBrewer, proposing a nice selection of color palettes. The viridis package, with many excellent palettes
#' @examples
#' fac2col(c("a", "b", "a", "c"))
fac2col <-
  function(x, col.pal = funky_color, na.col = "grey", seed = NULL) {
    x <- factor(x)
    lev <- levels(x)
    nlev <- length(lev)
    if (!is.null(seed)) {
      set.seed(seed)
      newseed <- round(runif(1, 1, 1e+09))
      on.exit(set.seed(newseed))
      col <- sample(col.pal(nlev))
    } else {
      col <- col.pal(nlev)
    }
    res <- rep(na.col, length(x))
    res[!is.na(x)] <- col[as.integer(x[!is.na(x)])]
    return(res)
  }
################################################################################

################################################################################
#' Adds transparency to a vector of colors
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' @param col a vector of colors
#' @param alpha (default 0.5) a numeric value between 0 and 1 representing the alpha coefficient; 0: total transparency; 1: no transparency.
#' @return a color vector
#' @author Thibaut Jombart in `adegenet` package
#' @export
#' @seealso The R package RColorBrewer, proposing a nice selection of color palettes. The viridis package, with many excellent palettes
#' @examples
#' transp("red")
#' transp(c("red", "blue"), alpha = 0.3)
transp <- function(col, alpha = 0.5) {
  res <-
    apply(col2rgb(col), 2, function(c) {
      rgb(c[1] / 255, c[2] / 255, c[3] / 255, alpha)
    })
  return(res)
}
################################################################################

################################################################################
#' Subsample a fastq file copying the n_seq first sequences in a given folder
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Useful to test a pipeline on small fastq files.
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
#' subsample_fastq(list_fastq_files(system.file("extdata", package = "MiscMetabar")),
#'   paste0(tempdir(), "/output_fastq"),
#'   n = 10
#' )
#' unlink(paste0(tempdir(), "/output_fastq"), recursive = TRUE)
#' }
subsample_fastq <- function(
  fastq_files,
  folder_output = "subsample",
  nb_seq = 1000
) {
  for (f in unlist(fastq_files)) {
    if (!dir.exists(folder_output)) {
      dir.create(folder_output)
    }
    writeLines(
      readLines(f, n = nb_seq * 4),
      con = paste0(
        folder_output,
        "/",
        basename(f)
      )
    )
  }
}

################################################################################

################################################################################
#' Test if cutadapt is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
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

is_cutadapt_installed <- function(
  args_before_cutadapt = "source ~/miniconda3/etc/profile.d/conda.sh && conda activate cutadaptenv && "
) {
  writeLines(
    paste0(args_before_cutadapt, " cutadapt -h"),
    paste0(tempdir(), "/script_cutadapt.sh")
  )
  cutadapt_error_or_not <- try(
    system(
      paste0("bash ", tempdir(), "/script_cutadapt.sh 2>&1"),
      intern = TRUE
    ),
    silent = TRUE
  )
  unlink(paste0(tempdir(), "/script_cutadapt.sh"))

  return(!inherits(cutadapt_error_or_not, "try-error"))
}

#' Test if falco is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
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
  return(
    !inherits(
      try(system(paste0(path, " 2>&1"), intern = TRUE), silent = TRUE),
      "try-error"
    )
  )
}

#' Test if swarm is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
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
  return(
    !inherits(
      try(system(paste0(path, " -h 2>&1"), intern = TRUE), silent = TRUE),
      "try-error"
    )
  )
}

#' Test if blastn is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: blastn) Path to blastn (NCBI BLAST+)
#' @export
#' @return A logical that say if blastn is installed.
#'
#' @examples
#' MiscMetabar::is_blastn_installed()
#' @author Adrien Taudière

is_blastn_installed <- function(path = "blastn") {
  return(
    !inherits(
      try(
        system2(path, "-version", stdout = TRUE, stderr = TRUE),
        silent = TRUE
      ),
      "try-error"
    )
  )
}

#' Test if MultiQC is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: multiqc) Path to MultiQC
#' @export
#' @return A logical that say if MultiQC is installed.
#'
#' @examples
#' MiscMetabar::is_multiqc_installed()
#' @author Adrien Taudière

is_multiqc_installed <- function(path = "multiqc") {
  return(
    !inherits(
      try(
        system2(path, "--version", stdout = TRUE, stderr = TRUE),
        silent = TRUE
      ),
      "try-error"
    )
  )
}

#' Test if vsearch is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' Useful for testthat and examples compilation for R CMD CHECK and
#'   test coverage
#'
#' @param path (default: [find_vsearch()]) Path to vsearch
#' @export
#' @return A logical that say if vsearch is install in
#'
#' @examples
#' MiscMetabar::is_vsearch_installed()
#' @author Adrien Taudière
#' @seealso [find_vsearch()], [install_vsearch()]

is_vsearch_installed <- function(path = find_vsearch()) {
  return(
    !inherits(
      try(
        system2(path, "--version", stdout = TRUE, stderr = TRUE),
        silent = TRUE
      ),
      "try-error"
    )
  )
}

#' Test if mumu is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
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
  return(
    !inherits(
      try(system(paste0(path, " 2>&1"), intern = TRUE), silent = TRUE),
      "try-error"
    )
  )
}


#' Test if krona is installed.
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
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
  return(inherits(
    try(system(paste0(path, " 2>&1"), intern = TRUE), silent = TRUE),
    "try-error"
  ))
}
################################################################################
