################################################################################
#' Make Krona files using [KronaTools](https://github.com/marbl/Krona/wiki).
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' Need the installation of kronatools on the computer ([installation instruction](https://github.com/marbl/Krona/wiki/Installing)).
#'
#' @param physeq (required) a \code{\link{phyloseq-class}} object.
#' @param file (required) the location of the html file to save
#' @param nb_seq (logical) If true, Krona set the distribution
#'   of sequences in the taxonomy. If False, Krona set the distribution of ASVs
#'   in the taxonomy.
#' @param ranks Number of the taxonomic ranks to plot
#'   (num of the column in `tax_table` slot of your `physeq` object).
#'   Default setting plot all the ranks (argument 'All').
#' @param add_unassigned_rank (int) Add unassigned for rank
#'   inferior to 'add_unassigned_rank' when necessary.
#' @param name A name for intermediary files, Useful to name
#'   your krona result files before merging using [merge_krona()]
#'
#' @examples
#' data("GlobalPatterns")
#' GA <- subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
#' # krona(GA, "Number.of.sequences.html")
#' # krona(GA, "Number.of.ASVs.html", nb_seq = F)
#' # merge_krona(c("Number.of.sequences.html", "Number.of.ASVs.html"))
#' @return A html file
#' @export
#' @seealso \code{\link{merge_krona}}
#' @author Adrien Taudière
#'

krona <-
  function(physeq,
           file = "krona.html",
           nb_seq = TRUE,
           ranks = "All",
           add_unassigned_rank = 0,
           name = NULL) {
    if (ranks[1] == "All") {
      ranks <- seq_along(physeq@tax_table[1, ])
    }

    df <- data.frame(unclass(physeq@tax_table[, ranks]))
    df$ASVs <- rownames(physeq@tax_table)

    if (is.null(name)) {
      if (nb_seq) {
        name <- "Number.of.sequences"
      } else {
        name <- "Number.of.ASVs"
      }
    }

    if (nb_seq) {
      df$nb_seq <- taxa_sums(physeq)
    } else {
      df$nb_otu <- rep(1, length(taxa_sums(physeq)))
    }

    df <- df[c(ncol(df), 2:ncol(df) - 1)]
    res <-
      lapply(split(df, seq_along(physeq@tax_table[, 1])), function(x) {
        as.vector(as.matrix(x))[!is.na(unlist(x))]
      })

    res <-
      lapply(res, function(x) {
        if (length(x) < add_unassigned_rank) {
          x <-
            c(x, "unassigned")[c(seq_len(x) - 1, length(x) + 1, length(x))]
        } else {
          x
        }
      })

    interm_txt <- paste(tempdir(), "/", name, ".html", sep = "")

    lapply(res,
      cat,
      "\n",
      file = interm_txt,
      append = TRUE,
      sep = "\t"
    )

    cmd <- paste("ktImportText ", interm_txt, " -o ", file, sep = "")
    system(command = cmd)
    system(command = paste("rm", interm_txt))
  }

###############################################################################
#' Merge Krona files using [KronaTools](https://github.com/marbl/Krona/wiki).
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' Need the installation of kronatools on the computer ([installation instruction](https://github.com/marbl/Krona/wiki/Installing)).
#' Function merge_krona allows merging multiple html files in one interactive krona file
#'
#' @param files (required) path to html files to merged
#' @param output path to the output file
#'
#' @examples
#' data("GlobalPatterns")
#' GA <- subset_taxa(GlobalPatterns, Phylum == "Acidobacteria")
#' # krona(GA, "Number.of.sequences.html")
#' # krona(GA, "Number.of.ASVs.html", nb_seq = F)
#' # merge_krona(c("Number.of.sequences.html", "Number.of.ASVs.html"))
#' @return A html file
#' @seealso \code{\link{krona}}
#' @export
#' @author Adrien Taudière
#'
merge_krona <- function(files = NULL, output = "mergeKrona.html") {
  cmd <-
    paste("ktImportKrona ",
      paste(files, collapse = " "),
      " -o ",
      output,
      sep = ""
    )

  system(command = cmd)
}
################################################################################
