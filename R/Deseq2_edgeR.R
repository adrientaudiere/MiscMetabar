################################################################################
#' Plot edgeR results for a phyloseq or a edgeR object.
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' @inheritParams clean_pq
#' @param contrast (required):This argument specifies what comparison
#'   to extract from the object to build a results table.
#'   See \code{\link[DESeq2]{results}} man page for more details.
#' @param pval (default: 0.05): the significance cutoff used for optimizing
#'   the independent filtering. If the adjusted p-value cutoff (FDR) will be a
#'   value other than 0.05, pval should be set to that value.
#' @param taxolev taxonomic level of interest
#' @param color_tax taxonomic level used for
#'   color assignation
#' @param verbose (logical): whether the function print some
#'   information during the computation
#' @param ... Additional arguments passed on to \code{\link[edgeR]{exactTest}}
#'   or \code{\link[ggplot2]{ggplot}}
#'
#' @export
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' data("GlobalPatterns", package = "phyloseq")
#' GP_archae <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' \donttest{
#' if (requireNamespace("edgeR")) {
#'   plot_edgeR_pq(GP_archae, c("SampleType", "Soil", "Feces"),
#'     color_tax = "Kingdom"
#'   )
#'
#'   plot_edgeR_pq(GP_archae, c("SampleType", "Soil", "Feces"),
#'     taxolev = "Class", color_tax = "Kingdom"
#'   )
#' }
#' }
#' @author Adrien Taudière
#'
#' @return A \code{\link{ggplot}}2 plot representing edgeR results
#'
#' @seealso \code{\link[edgeR]{exactTest}}
#' @seealso \code{\link{plot_deseq2_pq}}

plot_edgeR_pq <-
  function(physeq,
           contrast = NULL,
           pval = 0.05,
           taxolev = "Genus",
           color_tax = "Phylum",
           verbose = TRUE,
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("data must be an object of class 'phyloseq'")
    }

    if (verbose) {
      message("Conversion to edgeR format")
    }
    data_edger <- phyloseq_to_edgeR(physeq, group = contrast[1])
    if (verbose) {
      message("Perform edgeR binary test")
    }
    et <-
      edgeR::exactTest(data_edger, pair = c(contrast[2], contrast[3]), ...)

    tt <-
      edgeR::topTags(
        et,
        n = nrow(et$table),
        adjust.method = "BH",
        sort.by = "PValue"
      )
    res <- tt@.Data[[1]]
    sigtab <- res[(res$FDR < pval), ]
    sigtab <- cbind(methods::as(sigtab, "data.frame"))

    sigtabgen <- subset(sigtab, !is.na(taxolev))

    d <-
      tapply(sigtabgen$logFC, sigtabgen[, color_tax], function(x) {
        max(x)
      })
    d <- sort(d, TRUE)
    sigtabgen$col_tax <-
      factor(as.character(sigtabgen[, color_tax]), levels = names(d))

    d <-
      tapply(sigtabgen$logFC, sigtabgen[, taxolev], function(x) {
        max(x)
      })
    d <- sort(d, TRUE)
    sigtabgen$tax <-
      factor(as.character(sigtabgen[, taxolev]), levels = names(d))

    p <-
      ggplot(sigtabgen, aes(x = tax, y = logFC, color = col_tax), ...) +
      geom_point(size = 6) +
      theme(axis.text.x = element_text(
        angle = -90,
        hjust = 0,
        vjust = 0.5
      )) +
      labs(
        title = paste(
          "Change in abundance for ",
          contrast[1],
          " (",
          contrast[2],
          " vs ",
          contrast[3],
          ")",
          sep = ""
        )
      )

    return(p)
  }
################################################################################

################################################################################
# Plot the result of a DESeq2 test
################################################################################
#' Plot DESeq2 results for a phyloseq or a DESeq2 object.
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' @param data (required) a \code{\link{phyloseq-class}} or a
#' \code{\link[DESeq2]{DESeqDataSet-class}} object.
#' @param tax_table Required if data is a
#' \code{\link[DESeq2]{DESeqDataSet-class}} object.
#' The taxonomic table used to find the \code{taxa} and \code{color_taxa}
#' arguments. If data is a \code{\link{phyloseq-class}} object, data@tax_table
#' is used.
#' @param contrast (required) contrast specifies what comparison to extract
#'   from the object to build a results table. See \code{\link[DESeq2]{results}}
#'   man page for more details.
#' @param pval (default: 0.05) the significance cutoff used for optimizing
#'   the independent filtering. If the adjusted p-value cutoff (FDR) will be a
#'   value other than 0.05, pval should be set to that value.
#' @param taxolev taxonomic level of interest
#' @param select_taxa Either the name of the taxa (in the form of [DESeq2::results()])
#'   or a logical vector (length of the results from [DESeq2::results()]) to select taxa
#'   to plot.
#' @param color_tax taxonomic level used for color or a
#'   color vector.
#' @param tax_depth Taxonomic depth to test for differential
#'   distribution among contrast. If Null the analysis is done at the OTU
#'   (i.e. Species) level. If not Null, data need to be a column name in
#'   the `tax_table` slot of the \code{\link{phyloseq-class}} object.
#' @param verbose whether the function print some information during
#'   the computation
#' @param jitter_width width for the jitter positioning
#' @param ... Additional arguments passed on to \code{\link[DESeq2]{DESeq}}
#'   or \code{\link[ggplot2]{ggplot}}
#'
#' @export
#' @details
#'   Please cite `DESeq2` package if you use chis function.
#' @examples
#' \donttest{
#'
#' data("GlobalPatterns", package = "phyloseq")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' GP <- subset_samples(GP, SampleType %in% c("Soil", "Skin"))
#' if (requireNamespace("DESeq2")) {
#'   res <- DESeq2::DESeq(phyloseq_to_deseq2(GP, ~SampleType),
#'     test = "Wald", fitType = "local"
#'   )
#'   plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"),
#'     tax_table = GP@tax_table, color_tax = "Kingdom"
#'   )
#'   plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"),
#'     tax_table = GP@tax_table, color_tax = "Kingdom",
#'     pval = 0.7
#'   )
#'   plot_deseq2_pq(res, c("SampleType", "Soil", "Skin"),
#'     tax_table = GP@tax_table, color_tax = "Class",
#'     select_taxa = c("522457", "271582")
#'   )
#' }
#' }
#' @author Adrien Taudière
#'
#' @return A \code{\link{ggplot}}2 plot representing DESeq2 results
#'
#' @seealso \code{\link[DESeq2]{DESeq}}
#' @seealso \code{\link[DESeq2]{results}}
#' @seealso \code{\link{plot_edgeR_pq}}

plot_deseq2_pq <-
  function(data,
           contrast = NULL,
           tax_table = NULL,
           pval = 0.05,
           taxolev = "Genus",
           select_taxa = NULL,
           color_tax = "Phylum",
           tax_depth = NULL,
           verbose = TRUE,
           jitter_width = 0.1,
           ...) {
    if (!inherits(data, "phyloseq")) {
      if (!inherits(data, "DESeqDataSet")) {
        stop("data must be an object of class 'phyloseq' or 'DESeqDataSet'")
      }
    } else {
      # Calculate new dataset given the Taxa depth if tax_depth is not null
      if (!is.null(tax_depth)) {
        data_tax <- data
        data_tax@otu_table <-
          otu_table(
            apply(data@otu_table, 2, function(x) {
              tapply(
                x,
                data@tax_table[, tax_depth], sum
              )
            }),
            taxa_are_rows = TRUE
          )
        data_tax@tax_table <-
          tax_table(apply(
            data@tax_table[, 1:match(
              tax_depth,
              colnames(data@tax_table)
            )],
            2, function(x) {
              tapply(
                x, data@tax_table[, tax_depth],
                function(xx) {
                  xx[1]
                }
              )
            }
          ))
        data_tax@refseq <- NULL
        data <- data_tax
        if (is.na(match(taxolev, colnames(data@tax_table)))) {
          taxolev <- tax_depth
        }
      }

      if (is.null(tax_table) && inherits(data, "phyloseq")) {
        tax_table <- data@tax_table
      }

      if (verbose) {
        message("Conversion to Deseq2 format.")
      }
      data_deseq2 <-
        phyloseq_to_deseq2(data, stats::as.formula(paste("~", contrast[1])))

      if (verbose) {
        message("Calculation of Deseq2 results.")
      }
      data <-
        DESeq2::DESeq(
          data_deseq2,
          test = "Wald",
          fitType = "parametric",
          quiet = !verbose,
          ...
        )
    }

    # Calcul deseq2 results
    res <- DESeq2::results(data, contrast = contrast)

    if (!is.null(select_taxa[1])) {
      res <- res[select_taxa, ]
    }

    d <- res[which(res$padj < pval), ]

    if (dim(d)[1] == 0) {
      message("None taxa present significant distribution pattern through
              contrast.")
      return("None taxa present significant distribution pattern through
             contrast.")
    }
    d <-
      cbind(methods::as(d, "data.frame"), methods::as(tax_table[rownames(d), ], "matrix"))

    # Compute colors
    are_colors <- function(x) {
      sapply(x, function(xx) {
        tryCatch(
          is.matrix(grDevices::col2rgb(xx)),
          error = function(e) {
            FALSE
          }
        )
      })
    }

    if (!sum(are_colors(color_tax)) > 0) {
      x <- tapply(d$log2FoldChange, d[, color_tax], function(x) {
        max(x)
      })
      x <- sort(x, TRUE)
      d$col_tax <-
        factor(as.character(d[, color_tax]), levels = names(x))
    } else {
      d$col_tax <- rep(color_tax, length = dim(d)[1])
    }

    # Compute log2FoldChange values
    x <- tapply(d$log2FoldChange, d[, taxolev], function(x) {
      max(x)
    })
    x <- sort(x, TRUE)
    d$tax <- factor(as.character(d[, taxolev]), levels = names(x))

    if (!sum(are_colors(color_tax)) > 0) {
      p <-
        ggplot(d, aes(x = tax, y = log2FoldChange, color = col_tax), ...) +
        geom_point(
          size = 6,
          position = position_jitter(width = jitter_width, height = 0)
        ) +
        theme(axis.text.x = element_text(
          angle = -90,
          hjust = 0,
          vjust = 0.5
        )) +
        labs(
          title = paste(
            "Change in abundance for ",
            contrast[1],
            " (top:",
            contrast[2],
            " vs down:",
            contrast[3],
            ")",
            sep = ""
          )
        )
    } else {
      p <-
        ggplot(d, aes(x = tax, y = log2FoldChange), ...) +
        geom_point(
          size = 6, color = d$col_tax,
          position = position_jitter(width = jitter_width, height = 0)
        ) +
        theme(axis.text.x = element_text(
          angle = -90,
          hjust = 0,
          vjust = 0.5
        )) +
        labs(
          title = paste(
            "Change in abundance for ",
            contrast[1],
            " (top:",
            contrast[2],
            " vs down:",
            contrast[3],
            ")",
            sep = ""
          )
        )
    }

    return(p)
  }
################################################################################


################################################################################
#' Convert phyloseq OTU count data into DGEList for edgeR package
#'
#' @inheritParams clean_pq
#'
#' @param group (required) A character vector or factor giving the experimental
#'  group/condition for each sample/library. Alternatively, you may provide
#'  the name of a sample variable. This name should be among the output of
#'  \code{sample_variables(physeq)}, in which case
#'  \code{get_variable(physeq, group)} would return either a character vector or
#'    factor.
#'  This is passed on to \code{\link[edgeR]{DGEList}},
#'  and you may find further details or examples in its documentation.
#'
#' @param method The label of the edgeR-implemented normalization
#'   to use.
#'  See \code{\link[edgeR]{calcNormFactors}} for supported options and details.
#'  The default option is \code{"RLE"}, which is a scaling factor method
#'  proposed by Anders and Huber (2010).
#'  At time of writing, the \link[edgeR]{edgeR} package supported
#'  the following options to the \code{method} argument:
#'
#'  \code{c("TMM", "RLE", "upperquartile", "none")}.
#'
#' @return A DGEList object. See [edgeR::estimateTagwiseDisp()] for more details.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#' @export
#'
phyloseq_to_edgeR <- function(physeq, group, method = "RLE", ...) {
  verify_pq(physeq)
  # Enforce orientation.
  if (!taxa_are_rows(physeq)) {
    otu_table(physeq) <- otu_table(t(as.matrix(unclass(physeq@otu_table))),
      taxa_are_rows = TRUE
    )
  }
  x <- methods::as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x <- x + 1
  # Check `group` argument
  if (identical(all.equal(length(group), 1), TRUE) &&
    nsamples(physeq) > 1) {
    # Assume that group was a sample variable name (must be categorical)
    group <- get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy <- tax_table(physeq, errorIfNULL = FALSE)
  if (!is.null(taxonomy)) {
    taxonomy <- data.frame(methods::as(taxonomy, "matrix"))
  }
  # Now turn into a DGEList
  y <- edgeR::DGEList(
    counts = x,
    group = group,
    genes = taxonomy,
    remove.zeros = TRUE,
    ...
  )
  # Calculate the normalization factors
  z <- edgeR::calcNormFactors(y, method = method)
  # Check for division by zero inside `calcNormFactors`
  if (!all(is.finite(z$samples$norm.factors))) {
    stop(
      "Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument"
    )
  }
  # Estimate dispersions
  return(edgeR::estimateTagwiseDisp(edgeR::estimateCommonDisp(z)))
}
################################################################################
