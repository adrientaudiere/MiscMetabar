################################################################################
#' Plot the result of a mt test [phyloseq::mt()]
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Graphical representation of mt test.
#' 
#' @param mt (required) Result of a mt test from the function [phyloseq::mt()].
#' @param alpha (default: 0.05) Choose the cut off p-value to plot taxa.
#' @param color_tax (default: "Class") A taxonomic level to color the points.
#' @param taxa (default: "Species") The taxonomic level you choose for x-positioning.
#' @author Adrien Taudière
#' @examples
#' \donttest{
#' # Filter samples that don't have Time
#' data_fungi_mini2 <- subset_samples(data_fungi_mini, !is.na(Time))
#' res <- mt(data_fungi_mini2, "Time", method = "fdr", test = "f", B = 300)
#' plot_mt(res)
#' plot_mt(res, taxa = "Genus", color_tax = "Order")
#' }
#' @return a \code{\link{ggplot}}2 plot of result of a mt test
#' @export
#' @seealso [phyloseq::mt()]

plot_mt <-
  function(mt = NULL,
           alpha = 0.05,
           color_tax = "Class",
           taxa = "Species") {
    d <- mt[mt$plower < alpha, ]
    d$tax_col <- factor(as.character(d[, color_tax]))
    d$tax_col[is.na(d$tax_col)] <- "unidentified"
    d$tax <- as.character(d[, taxa])
    d$tax[is.na(d$tax)] <- "unidentified"
    d$tax <-
      factor(d$tax,
        levels =
          unique(factor(as.character(d[, taxa]))
          [rev(order(d$teststat))])
      )

    p <-
      ggplot(d, aes(x = tax, y = teststat, color = tax_col)) +
      geom_point(size = 6) +
      theme(axis.text.x = element_text(
        angle = -90,
        hjust = 0,
        vjust = 0.5
      ))
    p
  }
################################################################################


################################################################################
#' Plot accumulation curves for \code{\link{phyloseq-class}} object
#' 
#' @description 
#' 
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#' 
#' Note that as most bioinformatic pipeline discard singleton, accumulation curves from metabarcoding
#' cannot be interpreted in the same way as with conventional biodiversity sampling techniques.
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot
#'    different lines
#' @param add_nb_seq (default: TRUE, logical)
#' Either plot accumulation curves using sequences or using samples
#' @param step (Integer) distance among points calculated to plot lines. A
#'  low value give better plot but is more time consuming.
#'  Only used if `add_nb_seq` = TRUE.
#' @param by.fact (default: FALSE, logical)
#' First merge the OTU table by factor to plot only one line by factor
#' @param ci_col Color vector for confidence interval.
#'   Only use if `add_nb_seq` = FALSE.
#'   If `add_nb_seq` = TRUE, you can use ggplot to modify the plot.
#' @param col Color vector for lines. Only use if `add_nb_seq` = FALSE.
#'   If `add_nb_seq` = TRUE, you can use ggplot to modify the plot.
#' @param lwd  (default: 3) thickness for lines. Only use if `add_nb_seq` = FALSE.
#' @param leg (default: TRUE, logical) Plot legend or not. Only use if `add_nb_seq` = FALSE.
#' @param print_sam_names (default: FALSE, logical) Print samples names or not?
#'    Only use if `add_nb_seq` = TRUE.
#' @param ci (default: 2, integer) Confidence interval value used to multiply the
#'   standard error to plot confidence interval
#' @param ... Additional arguments passed on to \code{\link{ggplot}}
#' if `add_nb_seq` = TRUE or to \code{\link{plot}} if `add_nb_seq` = FALSE
#'
#' @examples
#' \donttest{
#' data("GlobalPatterns", package = "phyloseq")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' GP <- rarefy_even_depth(subset_samples_pq(GP, sample_sums(GP) > 3000))
#' p <- accu_plot(GP, "SampleType", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
#' p <- accu_plot(GP, "SampleType", add_nb_seq = TRUE, step = 10)
#'
#' p + theme(legend.position = "none")
#'
#' p + xlim(c(0, 400))
#' }
#' @return A \code{\link{ggplot}}2 plot representing the richness
#' accumulation plot if add_nb_seq = TRUE, else, if add_nb_seq = FALSE
#' return a base plot.
#'
#' @export
#' @author Adrien Taudière
#' @seealso \code{\link[vegan]{specaccum}} [accu_samp_threshold()]
accu_plot <-
  function(physeq,
           fact = NULL,
           add_nb_seq = TRUE,
           step = NULL,
           by.fact = FALSE,
           ci_col = NULL,
           col = NULL,
           lwd = 3,
           leg = TRUE,
           print_sam_names = FALSE,
           ci = 2,
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be a phyloseq object")
    }

    if (!taxa_are_rows(physeq)) {
      physeq@otu_table <-
        otu_table(t(physeq@otu_table), taxa_are_rows = TRUE)
    }

    if (!add_nb_seq) {
      factor_interm <-
        eval(parse(text = paste("physeq@sam_data$", fact, sep = "")))
      factor_interm <- as.factor(factor_interm)

      physeq_accu <- as.matrix(t(physeq@otu_table))
      physeq_accu[physeq_accu > 0] <- 1
      accu_all <- vegan::specaccum(physeq_accu)

      accu <- list()
      for (i in seq_along(levels(factor_interm))) {
        accu[[i]] <-
          vegan::specaccum(physeq_accu[factor_interm ==
            levels(factor_interm)[i], ])
      }

      if (is.null(col)) {
        col <- funky_color(nlevels(factor_interm) + 1)
      }
      if (is.null(ci_col)) {
        transp <- function(col, alpha = 0.5) {
          res <-
            apply(grDevices::col2rgb(col), 2, function(c) {
              grDevices::rgb(c[1] / 255, c[2] / 255, c[3] / 255, alpha)
            })
          return(res)
        }
        ci_col <-
          transp(funky_color(nlevels(factor_interm) + 1), 0.3)
      }

      plot(accu_all,
        # ci_type = "poly",
        # ci_col = ci_col[1],
        col = col[1],
        lwd = lwd,
        # ci_lty = 0,
        xlab = "Sample",
        ...
      )

      for (i in seq_along(levels(factor_interm))) {
        graphics::lines(accu[[i]],
          col = col[i + 1],
          lwd = lwd
        )
      }
      if (leg) {
        graphics::legend(
          "bottomright",
          c("all", levels(factor_interm)),
          col = col,
          lty = 1,
          lwd = 3
        )
      }
    }

    if (add_nb_seq) {
      fact_interm <-
        as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))

      if (!by.fact) {
        x <- t(physeq@otu_table)
      } else {
        x <- apply(physeq@otu_table, 1, function(x) {
          tapply(x, fact_interm, sum)
        })
      }

      tot <- rowSums(x)
      nr <- nrow(x)

      if (is.null(step)) {
        step <- round(max(tot) / 30, 0)
      }

      n_max <- seq(1, max(tot), by = step)
      out <- lapply(seq_len(nr), function(i) {
        res <- list()
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) {
          n <- c(n, tot[i])
        }
        res_interm <-
          vegan::rarefy(as.matrix(unclass(x[i, ])), n, se = TRUE)
        res <-
          cbind(as.matrix(res_interm)[1, ], as.matrix(res_interm)[2, ])
        return(res)
      })

      names(out) <- names(tot)

      df <- plyr::ldply(out, data.frame)

      cond <- vector(mode = "logical")
      for (i in seq_along(levels(as.factor(df$.id)))) {
        cond <- c(cond, 1:table(df$.id)[i])
      }

      df$x <- n_max[cond]

      if (!by.fact) {
        df$fact <-
          as.factor(unlist(unclass(physeq@sam_data
          [match(df$.id, sample_names(physeq)), fact])
          [fact]))
      } else {
        df$fact <- df$.id
      }

      df$ymin <- df$X1 - df$X2 * ci
      df$ymin[is.na(df$ymin)] <- df$X1[is.na(df$ymin)]
      df$ymax <- df$X1 + df$X2 * ci
      df$ymax[is.na(df$ymax)] <- df$X1[is.na(df$ymax)]
      dff <- data.frame(matrix(nrow = length(tot)))
      dff$xlab <- tapply(df$x, df$.id, max)
      dff$xlab <- dff$xlab + max(dff$xlab, na.rm = TRUE) / 20
      dff$ylab <- tapply(df$X1, df$.id, max)
      dff$.id <- names(dff$ylab)
      p <-
        ggplot(data = df, aes(
          x = x,
          y = X1,
          group = .id,
          col = fact
        )) +
        geom_ribbon(aes(
          ymin = ymin,
          ymax = ymax,
          col = NULL,
          fill = fact
        ), alpha = 0.2) +
        geom_line() +
        xlab("Number of sequences") +
        ylab("Number of OTUs (with standard error)")

      if (print_sam_names) {
        p + geom_text(data = dff, aes(
          x = xlab,
          y = ylab,
          label = .id,
          col = NULL
        ))
      } else {
        p
      }
      return(p)
    }
  }
################################################################################


################################################################################
#' Plot accumulation curves with balanced modality and depth rarefaction
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   This function (i) rarefy (equalize) the number of samples per modality of a
#'   factor and (ii) rarefy the number of sequences per sample (depth). The
#'   seed is set to 1:nperm. Thus, with exacly the same parameter, including
#'   nperm values, results must be identical.
#'
#' @inheritParams clean_pq
#' @param fact (required) The variable to rarefy. Must be present in
#'   the `sam_data` slot of the physeq object.
#' @param nperm (int) The number of permutations to perform.
#' @param step 	(int) distance among points calculated to plot lines.
#' A low value give better plot but is more time consuming.
#' @param by.fact (logical, default TRUE)
#' First merge the OTU table by factor to plot only one line by factor
#' @param progress_bar (logical, default TRUE) Do we print progress during
#'   the calculation?
#' @param quantile_prob (float, `[0:1]`) the value to compute the quantile.
#'   Minimum quantile is compute using 1-quantile_prob.
#' @param rarefy_by_sample_before_merging (logical, default TRUE):
#'    rarefy_by_sample_before_merging = FALSE is buggy for the moment.Please
#'    only use rarefy_by_sample_before_merging = TRUE
#' @param sample.size (int) A single integer value equal to the number of
#'   reads being simulated, also known as the depth. See
#'   [phyloseq::rarefy_even_depth()].
#' @param verbose (logical). If TRUE, print additional informations.
#' @param ... Other params for be passed on to [accu_plot()] function
#'
#' @export
#' @author Adrien Taudière
#' @seealso [accu_plot()], [rarefy_sample_count_by_modality()], [phyloseq::rarefy_even_depth()]
#'
#' @return A ggplot2 plot representing the richness accumulation plot
#' @examples
#' \donttest{
#' data_fungi_woNA4Time <-
#'   subset_samples(data_fungi, !is.na(Time))
#' data_fungi_woNA4Time@sam_data$Time <- paste0("time-", data_fungi_woNA4Time@sam_data$Time)
#' accu_plot_balanced_modality(data_fungi_woNA4Time, "Time", nperm = 3)
#'
#' data_fungi_woNA4Height <-
#'   subset_samples(data_fungi, !is.na(Height))
#' accu_plot_balanced_modality(data_fungi_woNA4Height, "Height", nperm = 3)
#' }
accu_plot_balanced_modality <- function(physeq,
                                        fact,
                                        nperm = 99,
                                        step = 2000,
                                        by.fact = TRUE,
                                        progress_bar = TRUE,
                                        quantile_prob = 0.975,
                                        rarefy_by_sample_before_merging = TRUE,
                                        sample.size = 1000,
                                        verbose = FALSE,
                                        ...) {
  if (rarefy_by_sample_before_merging) {
    p_for_dim <- accu_plot(
      rarefy_sample_count_by_modality(
        rarefy_even_depth(
          physeq,
          rngseed = 1,
          sample.size = sample.size,
          verbose = verbose
        ),
        fact,
        rngseed = 1
      ),
      fact = fact,
      step = step,
      by.fact = by.fact
    )$data
  } else {
    p_for_dim <- accu_plot(physeq,
      fact = fact,
      step = step,
      by.fact = by.fact
    )$data
    dim_for_plist <-
      max(tapply(sample_sums(physeq), physeq@sam_data[[fact]], sum))
  }

  dim_for_plist <- dim(p_for_dim)
  plist <- array(dim = c(dim_for_plist[1], 5, nperm))

  if (progress_bar) {
    pb <- txtProgressBar(
      min = 0,
      max = nperm,
      style = 3,
      width = 50,
      char = "="
    )
  }
  for (i in 1:nperm) {
    if (rarefy_by_sample_before_merging) {
      plist[, , i] <-
        as.matrix(suppressWarnings(suppressMessages(
          accu_plot(
            rarefy_sample_count_by_modality(
              rarefy_even_depth(
                physeq,
                rngseed = i,
                sample.size =
                  sample.size,
                verbose = verbose
              ),
              fact,
              rngseed = i
            ),
            fact = fact,
            step = step,
            by.fact = by.fact,
            ...
          )
        ))$data[, c(2:4, 6, 7)])
    } else {
      res_interm <-
        as.matrix(suppressWarnings(suppressMessages(
          accu_plot(
            rarefy_sample_count_by_modality(physeq,
              fact,
              rngseed = i,
              verbose = verbose
            ),
            fact = fact,
            step = step,
            by.fact = by.fact,
            ...
          )
        ))$data[, c(2:4, 6, 7)])
      plist[1:nrow(res_interm), , i] <- res_interm
    }
    if (progress_bar) {
      setTxtProgressBar(pb, i)
    }
  }

  res_mean <- data.frame(apply(plist, 1:2, mean, na.rm = TRUE))
  colnames(res_mean) <- colnames(p_for_dim[, c(2:4, 6, 7)])
  res_mean$fact <- p_for_dim$.id

  res_mean$X1_lim1 <- apply(plist, 1:2, function(x) {
    quantile(x, probs = quantile_prob, na.rm = TRUE)
  })[, 1]
  res_mean$X1_lim2 <- apply(plist, 1:2, function(x) {
    quantile(x, probs = 1 - quantile_prob, na.rm = TRUE)
  })[, 1]

  res_mean$factor <- p_for_dim$.id

  res_mean <- res_mean %>%
    dplyr::filter(!is.na(X1)) %>%
    arrange(X1)

  p <- ggplot(res_mean, aes(x = x, y = X1, color = factor)) +
    geom_line(linewidth = 1.5) +
    geom_ribbon(
      aes(
        ymin = X1_lim1,
        ymax = X1_lim2,
        fill = factor
      ),
      alpha = 0.2,
      linetype = 2,
      linewidth = 0.2
    )
  return(p)
}
################################################################################


################################################################################
#' Compute the number of sequence to obtain a given proportion of ASV in
#'  accumulation curves
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Note that as most bioinformatic pipeline discard singleton, accumulation curves from metabarcoding
#' cannot be interpreted in the same way as with conventional biodiversity sampling techniques.
#'
#' @param res_accuplot the result of the function accu_plot()
#' @param threshold the proportion of ASV to obtain in each samples
#'
#' @return a value for each sample of the number of sequences needed
#'   to obtain `threshold` proportion of the ASV
#'
#' @examples
#' \donttest{
#' data("GlobalPatterns", package = "phyloseq")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' GP <- rarefy_even_depth(subset_samples_pq(GP, sample_sums(GP) > 3000))
#' p <- accu_plot(GP, "SampleType", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
#'
#' val_threshold <- accu_samp_threshold(p)
#'
#' summary(val_threshold)
#'
#' # Plot the number of sequences needed to accumulate 0.95% of ASV in 50%, 75%
#' # and 100% of samples
#' p + geom_vline(xintercept = quantile(val_threshold, probs = c(0.50, 0.75, 1)))
#' }
#' @export
#' @author Adrien Taudière
#' @seealso [accu_plot()]
accu_samp_threshold <- function(res_accuplot, threshold = 0.95) {
  res <- list()
  for (id in unique(res_accuplot$data$.id)) {
    data <- res_accuplot$data %>% dplyr::filter(.id == id)
    proportion <- data$X1 / max(data$X1)
    res[[id]] <- data$x[proportion > threshold][1]
  }
  return(unlist(res))
}


################################################################################



################################################################################
#' Plot OTU circle for \code{\link{phyloseq-class}} object
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#' 
#' Graphical representation of distribution of taxa across a factor.
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor to cluster samples by modalities.
#'        Need to be in \code{physeq@sam_data}.
#' @param taxa (default: 'Order') Name of the taxonomic rank of interest
#' @param nproc (default 1)
#'   Set to number of cpus/processors to use for parallelization
#' @param add_nb_seq (default: TRUE) Represent the number of sequences or the
#'    number of OTUs (add_nb_seq = FALSE)
#' @param rarefy (logical) Does each samples modalities need to be rarefy in
#'               order to compare them with the same amount of sequences?
#' @param min_prop_tax (default: 0.01) The minimum proportion for taxa to be
#'                     plotted
#' @param min_prop_mod (default: 0.1) The minimum proportion for modalities
#'                     to be plotted
#' @param gap_degree Gap between two neighbour sectors.
#'                    It can be a single value or a vector. If it is a vector,
#'                     the first value corresponds to the gap after the first
#'                     sector.
#' @param start_degree The starting degree from which the circle begins to
#'   draw. Note this degree is measured in the standard polar coordinate
#'   which means it is always reverse-clockwise.
#' @param row_col Color vector for row
#' @param grid_col Grid colors which correspond to sectors. The length of the
#'    vector should be either 1 or the number of sectors.
#'    It's preferred that grid_col is a named vector of which names
#'    correspond to sectors. If it is not a named vector, the
#'    order of grid_col corresponds to order of sectors.
#' @param log10trans (logical) Should sequence be log10 transformed
#'                   (more precisely by log10(1+x))?
#' @param ... Additional arguments passed on to
#'   \code{\link[circlize]{chordDiagram}} or \code{\link[circlize]{circos.par}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("pbapply")) {
#'   data("GlobalPatterns", package = "phyloseq")
#'   GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#'   circle_pq(GP, "SampleType")
#'   circle_pq(GP, "SampleType", add_nb_seq = FALSE)
#'   circle_pq(GP, "SampleType", taxa = "Class")
#' }
#' }
#' @author Adrien Taudière
#'
#' @return A \code{\link[circlize]{chordDiagram}} plot representing the
#'   distribution of OTUs or sequences in the different modalities of the factor
#'   fact
#'
#' @export
#' @seealso \code{\link[circlize]{chordDiagram}}
#' @seealso \code{\link[circlize]{circos.par}}

circle_pq <-
  function(physeq = NULL,
           fact = NULL,
           taxa = "Order",
           nproc = 1,
           add_nb_seq = TRUE,
           rarefy = FALSE,
           min_prop_tax = 0.01,
           min_prop_mod = 0.1,
           gap_degree = NULL,
           start_degree = NULL,
           row_col = NULL,
           grid_col = NULL,
           log10trans = FALSE,
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }

    if (!physeq@otu_table@taxa_are_rows) {
      otu_tab <- t(physeq@otu_table)
    } else {
      otu_tab <- physeq@otu_table
    }

    if (!add_nb_seq) {
      otu_tab[otu_tab > 0] <- 1
    }

    taxcol <- match(taxa, colnames(physeq@tax_table))
    if (is.na(taxcol)) {
      stop("The taxa argument do not match any taxa rank in physeq@tax_table")
    }

    taxsamp <- match(fact, colnames(physeq@sam_data))
    if (is.na(taxsamp)) {
      stop("The samples argument do not match any sample attributes
           in physeq@sam_data")
    }

    otu_table_tax <-
      pbapply::pbapply(otu_tab, 2, function(x) {
        pbapply::pbtapply(
          x, physeq@tax_table[, taxcol],
          function(xx) {
            sum(xx, na.rm = TRUE)
          }
        )
      }, cl = nproc)
    otu_table_ech <-
      pbapply::pbapply(otu_table_tax, 1, function(x) {
        pbapply::pbtapply(
          x, physeq@sam_data[, taxsamp],
          function(xx) {
            sum(xx, na.rm = TRUE)
          }
        )
      }, cl = nproc)
    if (rarefy) {
      otu_table_ech_interm <-
        vegan::rrarefy(otu_table_ech, min(rowSums(otu_table_ech)))
      message(
        paste(
          "Rarefaction by modalities deletes ",
          sum(otu_table_ech) - sum(otu_table_ech_interm),
          " (",
          round(
            100 * (
              sum(otu_table_ech) - sum(otu_table_ech_interm)
            ) / sum(otu_table_ech),
            2
          ),
          "%) sequences.",
          sep = ""
        )
      )
      otu_table_ech <- otu_table_ech_interm
    }

    otu_table_ech <- otu_table_ech[, colSums(otu_table_ech) > 0]

    # Keep only taxa and modalities with a sufficient proportion (min_prop_tax,
    # min_prop_mod) to plot
    o_t_e_interm <-
      otu_table_ech[
        (rowSums(otu_table_ech) / sum(otu_table_ech)) >
          min_prop_mod,
        (colSums(otu_table_ech) / sum(otu_table_ech)) >
          min_prop_tax
      ]
    if (nrow(o_t_e_interm) != nrow(otu_table_ech)) {
      message(
        paste(
          "Only ",
          nrow(o_t_e_interm),
          " modalities are plot (",
          round(100 *
            nrow(o_t_e_interm) /
            nrow(otu_table_ech), 2),
          "%). Use 'min_prop_mod' to plot more samples.",
          sep = ""
        )
      )
    }

    if (ncol(o_t_e_interm) != ncol(otu_table_ech)) {
      message(
        paste(
          "Only ",
          ncol(o_t_e_interm),
          " taxa are plot (",
          round(
            100 * ncol(o_t_e_interm) / ncol(otu_table_ech),
            2
          ),
          "%). Use 'min_prop_tax' to plot more taxa",
          sep = ""
        )
      )
    }
    otu_table_ech <- o_t_e_interm

    if (log10trans) {
      otu_table_ech <- apply(otu_table_ech, 2, function(x) {
        log10(1 + x)
      })
    }


    if (is.null(gap_degree)) {
      col2keep <- rep(1, ncol(otu_table_ech) - 1)
      row2keep <- rep(1, nrow(otu_table_ech) - 1)
      gap_degree <- c(row2keep, 10, col2keep, 10)
    }
    if (is.null(start_degree)) {
      start_degree <- 170
    }

    if (is.null(grid_col)) {
      grid_col <-
        c(funky_color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
    }

    if (is.null(row_col)) {
      row_col <-
        c(funky_color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
    }

    circlize::circos.par(
      gap.degree = gap_degree, start.degree = start_degree,
      ...
    )
    circlize::chordDiagram(otu_table_ech,
      row.col = row_col,
      grid.col = grid_col, ...
    )
    circlize::circos.clear()
  }
################################################################################

################################################################################
#' Sankey plot of \code{\link{phyloseq-class}} object
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Graphical representation of distribution of taxa across Taxonomy and (optionnaly a factor).
#' 
#' @inheritParams clean_pq
#' @param fact Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}.
#' @param taxa a vector of taxonomic rank to plot
#' @param add_nb_seq Represent the number of sequences or the
#'   number of OTUs (add_nb_seq = FALSE). Note that plotting the number of
#'   sequences is slower.
#' @param min_prop_tax (default: 0) The minimum proportion for taxa to be
#'  plotted. EXPERIMENTAL. For the moment each links below the min.prop.
#'  tax is discard from the sankey network resulting in sometimes weird plot.
#' @param tax2remove  a vector of taxonomic groups to remove from the analysis
#'   (e.g. \code{c('Incertae sedis', 'unidentified')})
#' @param units  character string describing physical units (if any) for Value
#' @param symbol2sub (default: c('\\.', '-')) vector of symbol to delete in
#'   the taxonomy
#' @param ... Additional arguments passed on to
#'   \code{\link[networkD3]{sankeyNetwork}}
#'
#' @examples
#' data("GlobalPatterns", package = "phyloseq")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == "Archaea")
#' if (requireNamespace("networkD3")) {
#'   sankey_pq(GP, fact = "SampleType")
#' }
#' \donttest{
#' if (requireNamespace("networkD3")) {
#'   sankey_pq(GP, taxa = 1:4, min_prop_tax = 0.01)
#'   sankey_pq(GP, taxa = 1:4, min_prop_tax = 0.01, add_nb_seq = TRUE)
#' }
#' }
#' @author Adrien Taudière
#'
#' @return A \code{\link[networkD3]{sankeyNetwork}} plot representing the
#'  taxonomic distribution of OTUs or sequences. If \code{fact} is set,
#'  represent the distribution of the last taxonomic level in the modalities
#'  of \code{fact}
#'
#' @export
#' @seealso \code{\link[networkD3]{sankeyNetwork}}, [ggaluv_pq()]

sankey_pq <-
  function(physeq = NULL,
           fact = NULL,
           taxa = 1:4,
           add_nb_seq = FALSE,
           min_prop_tax = 0,
           tax2remove = NULL,
           units = NULL,
           symbol2sub = c("\\.", "-"),
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }

    if (!physeq@otu_table@taxa_are_rows) {
      otu_tab <- t(physeq@otu_table)
    } else {
      otu_tab <- physeq@otu_table
    }

    if (!add_nb_seq) {
      otu_tab[otu_tab > 0] <- 1
      mat_interm <- matrix()
      mat <- matrix(ncol = 3)
      colnames(mat) <- c("Var1", "Var2", "value")
      for (i in 1:(length(taxa) - 1)) {
        res_interm <-
          table(physeq@tax_table[, taxa[i]], physeq@tax_table[, taxa[i + 1]])
        mat_interm <- reshape2::melt(res_interm)
        mat_interm <- mat_interm[mat_interm[, 3] > 0, ]
        mat <- rbind(mat, mat_interm)
      }
    } else if (add_nb_seq) {
      mat_interm <- matrix()
      mat <- matrix(ncol = 3)
      colnames(mat) <- c("Var1", "Var2", "value")
      tax_table_interm <-
        physeq@tax_table[rep(seq(1, ntaxa(physeq)),
          times = taxa_sums(physeq)
        )]

      for (i in 1:(length(taxa) - 1)) {
        res_interm <-
          table(tax_table_interm[, taxa[i]], tax_table_interm[, taxa[i + 1]])
        mat_interm <- reshape2::melt(res_interm)
        mat_interm <- mat_interm[mat_interm[, 3] > 0, ]
        mat <- rbind(mat, mat_interm)
      }
    }

    if (!is.null(fact)) {
      net_matrix2links <- function(m = NULL) {
        res <- matrix(ncol = 3)
        for (i in seq_len(dim(m)[1])) {
          for (j in seq_len(dim(m)[2])) {
            if (m[i, j] > 0) {
              res <- rbind(res, c(rownames(m)[i], colnames(m)[j], m[i, j]))
            }
          }
        }
        return(res)
      }

      mat_interm <-
        apply(otu_tab, 1, function(x) {
          tapply(
            x, physeq@sam_data[, fact],
            sum
          )
        })

      if (!add_nb_seq) {
        mat_interm <-
          apply(mat_interm, 1, function(x) {
            tapply(x, physeq@tax_table[
              ,
              taxa[length(taxa)]
            ], function(x) {
              sum(x > 0)
            })
          })
      } else if (add_nb_seq) {
        mat_interm <-
          apply(mat_interm, 1, function(x) {
            tapply(x, physeq@tax_table[
              ,
              taxa[length(taxa)]
            ], sum)
          })
      }

      samp_links <- net_matrix2links(mat_interm)
      samp_links[, 2] <- toupper(samp_links[, 2])
      colnames(samp_links) <- colnames(mat)
      mat <- rbind(mat, samp_links)
    }

    mat <- as.data.frame(mat[rowSums(is.na(mat)) == 0, ])
    mat[, 3] <- as.numeric(as.vector(mat[, 3]))
    mat <- mat[rowSums(is.na(mat)) == 0, ]


    if (!is.null(tax2remove)) {
      mat <- mat[!mat[, 1] %in% tax2remove, ]
      mat <- mat[!mat[, 2] %in% tax2remove, ]
    }

    if (min_prop_tax != 0) {
      min_nb_tax <- min_prop_tax * sum(mat[, 3]) / length(taxa)
      mat <- mat[mat[, 3] >= min_nb_tax, ]
    }

    for (i in seq_len(length(symbol2sub))) {
      mat <- apply(mat, 2, function(x) {
        gsub(symbol2sub[i], "", x)
      })
    }

    tax_sank <- list()
    names_nodes <-
      unique(c(as.vector(mat[, 1]), as.vector(mat[, 2])))
    names_nodes <- names_nodes[!is.na(names_nodes)]
    tax_sank$nodes <-
      data.frame((seq_len(length(names_nodes))) - 1, names_nodes)
    names(tax_sank$nodes) <- c("code", "name")
    mat2 <- mat
    for (i in seq_len(nrow(tax_sank$nodes))) {
      mat2[, 1] <-
        gsub(
          paste("\\<", tax_sank$nodes[i, 2], "\\>", sep = ""),
          tax_sank$nodes[
            i,
            1
          ],
          mat2[, 1]
        )
      mat2[, 2] <-
        gsub(
          paste("\\<", tax_sank$nodes[i, 2], "\\>", sep = ""),
          tax_sank$nodes[
            i,
            1
          ],
          mat2[, 2]
        )
    }

    tax_sank$links <- apply(mat2, 2, as.numeric)
    tax_sank$links <-
      data.frame(tax_sank$links[rowSums(is.na(tax_sank$links)) == 0, ])
    tax_sank$nodes <-
      as.data.frame(as.character(tax_sank$nodes[, 2]))
    names(tax_sank$nodes) <- "name"
    names(tax_sank$links) <- c("source", "target", "value")
    if (is.null(units)) {
      if (!add_nb_seq) {
        units <- "OTUs"
      } else if (add_nb_seq) {
        units <- "Sequences"
      }
    }
    networkD3::sankeyNetwork(
      Links = tax_sank$links,
      Nodes = tax_sank$nodes,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      units = units,
      ...
    )
  }
################################################################################

################################################################################
#' Venn diagram of \code{\link{phyloseq-class}} object
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#' 
#' Graphical representation of distribution of taxa across combined modality of a factor.
#'
#' @inheritParams clean_pq
#' @param fact (required): Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}.
#' @param min_nb_seq (default: 0)): minimum number of sequences by OTUs by
#'  samples to take into count this OTUs in this sample. For example,
#'  if min_nb_seq=2,each value of 2 or less in the OTU table
#'  will be change into 0 for the analysis
#' @param print_values (logical) Print (or not) the table of number of OTUs
#' for each combination.
#' If print_values is TRUE the object is not a ggplot object.
#' Please use print_values = FALSE if you want to add ggplot function
#' (cf example).
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' if (requireNamespace("venneuler")) {
#'   data("enterotype")
#'   venn_pq(enterotype, fact = "SeqTech")
#' }
#' \donttest{
#' if (requireNamespace("venneuler")) {
#'   venn_pq(enterotype, fact = "ClinicalStatus")
#'   venn_pq(enterotype, fact = "Nationality", print_values = FALSE)
#'   venn_pq(enterotype, fact = "ClinicalStatus", print_values = FALSE) +
#'     scale_fill_hue()
#'   venn_pq(enterotype, fact = "ClinicalStatus", print_values = FALSE) +
#'     scale_fill_hue()
#' }
#' }
#' @return A \code{\link{ggplot}}2 plot representing Venn diagram of
#' modalities of the argument \code{factor}
#'
#' @export
#' @author Adrien Taudière
#' @seealso \code{\link[venneuler]{venneuler}}

venn_pq <-
  function(physeq,
           fact,
           min_nb_seq = 0,
           print_values = TRUE) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }

    moda <-
      as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))
    if (length(moda) != dim(physeq@otu_table)[1]) {
      data_venn <-
        t(apply(physeq@otu_table, 1, function(x) {
          by(x, moda, max)
        }))
    } else if (length(moda) != dim(physeq@otu_table)[2]) {
      data_venn <-
        t(apply(t(physeq@otu_table), 1, function(x) {
          by(x, moda, max)
        }))
    } else {
      stop("The factor length and the number of samples must be identical")
    }
    combinations <- data_venn > min_nb_seq

    e <- new.env(TRUE, emptyenv())
    cn <- colnames(combinations)
    for (i in seq.int(dim(combinations)[1])) {
      if (any(combinations[i, ])) {
        ec <- paste(cn[combinations[i, ]], collapse = "&")
        e[[ec]] <- if (is.null(e[[ec]])) {
          1L
        } else {
          (e[[ec]] + 1L)
        }
      }
    }

    en <- ls(e, all.names = TRUE)
    weights <- as.numeric(unlist(lapply(en, get, e)))
    combinations <- as.character(en)

    table_value <-
      data.frame(
        combinations = as.character(combinations),
        weights = as.double(weights)
      )

    venn <- venneuler::venneuler(data_venn > min_nb_seq)
    venn_res <-
      data.frame(
        x = venn$centers[, 1],
        y = venn$centers[, 2],
        radius = venn$diameters / 2
      )

    nmod <- nrow(venn_res)
    x1 <- list()
    for (i in seq(1, nmod)) {
      x1[[i]] <- grep(rownames(venn_res)[i], table_value$combinations)
    }

    for (i in seq_len(nrow(table_value))) {
      table_value$x[i] <-
        mean(venn$centers[, "x"][unlist(lapply(
          x1,
          function(x) {
            sum(x %in% i) > 0
          }
        ))])
      table_value$y[i] <-
        mean(venn$centers[, "y"][unlist(lapply(
          x1,
          function(x) {
            sum(x %in% i) > 0
          }
        ))])
    }

    df <- venn_res
    df$xlab <- df$x + (df$x - mean(df$x))
    df$ylab <- df$y + (df$y - mean(df$y))

    circularise <- function(d, n = 360) {
      angle <- seq(-pi, pi, length = n)
      make_circle <- function(x, y, r, modality) {
        data.frame(
          x = x + r * cos(angle),
          y = y + r * sin(angle),
          modality
        )
      }
      lmat <- mapply(
        make_circle,
        modality = rownames(d),
        x = d[, 1],
        y = d[, 2],
        r = d[, 3],
        SIMPLIFY = FALSE
      )
      do.call(rbind, lmat)
    }

    circles <- circularise(df)

    p <-
      ggplot() +
      geom_polygon(
        data = circles,
        aes(x, y, group = modality, fill = modality),
        alpha = 0.5
      ) +
      theme_void()

    if (print_values) {
      g_legend <- function(agplot) {
        tmp <- ggplot_gtable(ggplot_build(agplot))
        leg <-
          which(vapply(tmp$grobs, function(x) {
            x$name
          }, character(1)) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }
      legend <- g_legend(p)


      grid::grid.newpage()
      vp1 <- grid::viewport(
        width = 0.75,
        height = 1,
        x = 0.375,
        y = .5
      )
      vpleg <-
        grid::viewport(
          width = 0.25,
          height = 0.5,
          x = 0.85,
          y = 0.75
        )
      subvp <- grid::viewport(
        width = 0.3,
        height = 0.3,
        x = 0.85,
        y = 0.25
      )
      print(p + theme(legend.position = "none"), vp = vp1)
      grid::upViewport(0)
      grid::pushViewport(vpleg)
      grid::grid.draw(legend)
      grid::upViewport(0)
      grid::pushViewport(subvp)
      grid::grid.draw(gridExtra::tableGrob(table_value[, c(1, 2)], rows = NULL))
    } else {
      return(p)
    }
  }
################################################################################



################################################################################
#' Venn diagram of \code{\link{phyloseq-class}} object using
#' `ggVennDiagram::ggVennDiagram` function
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Note that you can use ggplot2 function to customize the plot
#' for ex. `+ scale_fill_distiller(palette = "BuPu", direction = 1)`
#' and `+ scale_x_continuous(expand = expansion(mult = 0.5))`. See
#' examples.
#'
#' @inheritParams clean_pq
#' @param fact (required): Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}.
#' @param min_nb_seq minimum number of sequences by OTUs by
#'   samples to take into count this OTUs in this sample. For example,
#'   if min_nb_seq=2,each value of 2 or less in the OTU table
#'   will not count in the venn diagram
#' @param taxonomic_rank Name (or number) of a taxonomic rank
#'   to count. If set to Null (the default) the number of OTUs is counted.
#' @param split_by Split into multiple plot using variable split_by.
#'   The name of a variable must be present in `sam_data` slot
#'   of the physeq object.
#' @param add_nb_samples (logical, default TRUE) Add the number of samples to
#'    levels names
#' @param add_nb_seq (logical, default FALSE) Add the number of sequences to
#'    levels names
#' @param rarefy_before_merging Rarefy each sample before merging by the
#'   modalities of args `fact`. Use `phyloseq::rarefy_even_depth()` function
#' @param rarefy_after_merging Rarefy each sample after merging by the
#'   modalities of args `fact`.
#' @param ... Other arguments for the `ggVennDiagram::ggVennDiagram` function
#'   for ex. `category.names`.
#' @return A \code{\link{ggplot}}2 plot representing Venn diagram of
#'   modalities of the argument \code{factor} or if split_by is set a list
#'   of plots.
#' @seealso [upset_pq()]
#' @examples
#' if (requireNamespace("ggVennDiagram")) {
#'   ggvenn_pq(data_fungi, fact = "Height")
#' }
#' \donttest{
#' if (requireNamespace("ggVennDiagram")) {
#'   ggvenn_pq(data_fungi, fact = "Height") +
#'     ggplot2::scale_fill_distiller(palette = "BuPu", direction = 1)
#'   pl <- ggvenn_pq(data_fungi, fact = "Height", split_by = "Time")
#'   for (i in seq_along(pl)) {
#'     p <- pl[[i]] +
#'       scale_fill_distiller(palette = "BuPu", direction = 1) +
#'       theme(plot.title = element_text(hjust = 0.5, size = 22))
#'     print(p)
#'   }
#'
#'   data_fungi2 <- subset_samples(data_fungi, data_fungi@sam_data$Tree_name == "A10-005" |
#'     data_fungi@sam_data$Height %in% c("Low", "High"))
#'   ggvenn_pq(data_fungi2, fact = "Height")
#'
#'   ggvenn_pq(data_fungi, fact = "Height", add_nb_seq = TRUE, set_size = 4)
#'   ggvenn_pq(data_fungi, fact = "Height", rarefy_before_merging = TRUE)
#'   ggvenn_pq(data_fungi, fact = "Height", rarefy_after_merging = TRUE) +
#'     scale_x_continuous(expand = expansion(mult = 0.5))
#' }
#' }
#' @export
#' @author Adrien Taudière


ggvenn_pq <- function(physeq = NULL,
                      fact = NULL,
                      min_nb_seq = 0,
                      taxonomic_rank = NULL,
                      split_by = NULL,
                      add_nb_samples = TRUE,
                      add_nb_seq = FALSE,
                      rarefy_before_merging = FALSE,
                      rarefy_after_merging = FALSE,
                      ...) {
  if (!is.factor(physeq@sam_data[[fact]])) {
    physeq@sam_data[[fact]] <- as.factor(physeq@sam_data[[fact]])
  }

  physeq <- taxa_as_columns(physeq)

  if (rarefy_before_merging) {
    physeq <- rarefy_even_depth(physeq)
    physeq <- clean_pq(physeq)
  }

  nb_samples <- table(physeq@sam_data[[fact]])

  if (rarefy_after_merging) {
    physeq <- merge_samples2(physeq, fact)
    physeq <- rarefy_even_depth(physeq)
    physeq <- clean_pq(physeq)
  }

  res <- list()
  nb_seq <- vector(mode = "integer")

  for (f in levels(physeq@sam_data[[fact]])) {
    newphyseq <- physeq
    new_DF <- newphyseq@sam_data[newphyseq@sam_data[[fact]] == f, ]
    sample_data(newphyseq) <- sample_data(new_DF)
    if (is.null(taxonomic_rank)) {
      res[[f]] <- colnames(newphyseq@otu_table[
        ,
        colSums(newphyseq@otu_table) > min_nb_seq
      ])
    } else {
      res[[f]] <-
        as.character(stats::na.exclude(unique(newphyseq@tax_table[
          colSums(newphyseq@otu_table) > min_nb_seq,
          taxonomic_rank
        ])))
    }
    nb_seq <-
      c(nb_seq, sum(physeq@otu_table[physeq@sam_data[[fact]] == f, ], na.rm = TRUE))
  }

  if (max(nb_seq) / min(nb_seq) > 2) {
    message(
      paste0(
        "Two modalities differ greatly (more than x2) in their number of sequences (",
        max(nb_seq),
        " vs ",
        min(nb_seq),
        "). You may be interested by the parameter rarefy_after_merging"
      )
    )
  }

  if (add_nb_samples) {
    names(res) <- paste0(names(res), "\n (", nb_samples, " sam.)")
  }

  if (add_nb_seq) {
    names(res) <- paste0(names(res), "\n (", nb_seq, " seq.)")
  }

  if (is.null(split_by)) {
    p <- ggVennDiagram::ggVennDiagram(res, ...)
  } else {
    p <- list()
    modalities <-
      as.factor(unlist(unclass(physeq@sam_data[[split_by]])))
    for (moda in levels(modalities)) {
      physeq_interm <-
        clean_pq(subset_samples_pq(physeq, modalities == moda),
          silent = TRUE
        )
      p[[moda]] <- ggvenn_pq(
        physeq_interm,
        fact = fact,
        min_nb_seq = 0,
        taxonomic_rank = NULL
      ) +
        ggtitle(moda)
    }
  }
  return(p)
}
################################################################################




################################################################################
#' Multiple plot function
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-stable-green" alt="lifecycle-stable"></a>
#'
#' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot
#' objects)
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#'
#' @param ... list of ggplot objects
#' @param plotlist list of ggplot objects
#' @param cols number of columns
#' @param layout A matrix specifying the layout.
#'   If present, 'cols' is ignored.
#' @return Nothing. Print the list of ggplot objects
#' @export

multiplot <-
  function(...,
           plotlist = NULL,
           cols = 1,
           layout = NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    num_plots <- length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(num_plots / cols)),
        ncol = cols,
        nrow = ceiling(num_plots / cols)
      )
    }

    if (num_plots == 1) {
      message(plots[[1]])
    } else {
      # Set up the page
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(
        nrow(layout),
        ncol(layout)
      )))

      # Make each plot, in the correct location
      for (i in seq(1, num_plots)) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <-
          as.data.frame(which(layout == i, arr.ind = TRUE))

        print(
          plots[[i]],
          vp = grid::viewport(
            layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col
          )
        )
      }
    }
  }
################################################################################

################################################################################
#' Graphical representation of hill number 0, 1 and 2 across a factor
#' 
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Hill numbers are the number of equiprobable species giving the same
#'   diversity value as the observed distribution. The Hill number 0
#'   correspond to Species richness), the Hill number 1 to
#'   the exponential of Shannon Index and the Hill number 2 to the inverse
#'   of Simpson Index)
#'
#' Note that (if correction_for_sample_size is TRUE, default behavior)
#'   this function use a sqrt of the read numbers in the linear
#'   model in order to correct for uneven sampling depth. This correction
#'   is only done before tuckey HSD plot and do not change the hill number
#'   computed.
#' 
#' @inheritParams clean_pq
#' @param fact (required): The variable to test. Must be present in
#'   the `sam_data` slot of the physeq object.
#' @param variable : Alias for factor. Kept only for backward compatibility.
#' @param hill_scales (a vector of integer) The list of q values to compute
#'   the hill number H^q. If Null, no hill number are computed. Default value
#'   compute the Hill number 0 (Species richness), the Hill number 1
#'   (exponential of Shannon Index) and the Hill number 2 (inverse of Simpson
#'   Index).
#' @param color_fac (optional): The variable to color the barplot. For ex.
#'   same as fact. Not very useful because ggplot2 plot colors can be
#'   change using `scale_color_XXX()` function.
#' @param letters (optional, default FALSE): If set to TRUE, the plot
#'   show letters based on p-values for comparison. Use the
#'   \code{\link[multcompView]{multcompLetters}} function from the package
#'   multcompLetters. BROKEN for the moment. Note that na values in The
#'   variable param need to be removed (see examples) to use letters.
#' @param add_points (logical, default FALSE): add jitter point on boxplot
#' @param add_info (logical, default TRUE) Do we add a subtitle with
#'   information about the number of samples per modality ?
#' @param one_plot (logical, default FALSE) If TRUE, return a unique
#'   plot with the four plot inside using the patchwork package.
#'   Note that if letters and one_plot are both TRUE, tuckey HSD results
#'   are discarded from the unique plot. In that case, use one_plot = FALSE
#'   to see the tuckey HSD results in the fourth plot of the resulting list.
#' @param kruskal_test (logical, default TRUE) Do we test for global effect of
#'   our factor on each hill scales values? When kruskal_test is TRUE, the
#'   resulting test value are add in each plot in subtitle (unless add_info is
#'   FALSE). Moreover, if at
#'   least one hill scales is not significantly link to fact (pval>0.05),
#'   a message is prompt saying that Tuckey HSD plot is not informative for
#'   those Hill scales and letters are not printed.
#' @param plot_with_tuckey (logical, default TRUE). If one_plot is set to
#'   TRUE and letters to FALSE, allow to discard the tuckey plot part with
#'   plot_with_tuckey = FALSE
#' @param correction_for_sample_size (logical, default TRUE) This function
#'   use a sqrt of the read numbers in the linear model in order to
#'   correct for uneven sampling depth in the Tuckey TEST. This params
#'   do not change value of Hill number but only the test associated
#'   values (including the pvalues). To rarefy samples, you may use the
#'   function [phyloseq::rarefy_even_depth()].
#' @param na_remove (logical, default TRUE) Do we remove samples with NA in
#'   the factor fact ? Note that na_remove is always TRUE when using
#'   letters = TRUE
#' @param vioplot (logical, default FALSE) Do we plot violin plot instead of
#'   boxplot ?
#' @return Either an unique ggplot2 object (if one_plot is TRUE) or
#'  a list of n+1 ggplot2 plot (with n the number of hill scale value).
#'  For example, with the default scale value:
#' - plot_Hill_0 : the boxplot of Hill number 0 (= species richness)
#'     against the variable
#' - plot_Hill_1 : the boxplot of Hill number 1 (= Shannon index)
#'      against the variable
#' - plot_Hill_2 : the boxplot of Hill number 2 (= Simpson index)
#'     against the variable
#' - plot_tuckey : plot the result of the Tuckey HSD test
#'
#' @export
#' @author Adrien Taudière
#' @examples
#'
#' p <- hill_pq(data_fungi_mini, "Height", hill_scales = 1:2)
#' p_h1 <- p[[1]] + theme(legend.position = "none")
#' p_h2 <- p[[2]] + theme(legend.position = "none")
#' multiplot(plotlist = list(p_h1, p_h2, p[[3]]), cols = 4)
#' \donttest{
#' if (requireNamespace("multcompView")) {
#'   p2 <- hill_pq(data_fungi, "Time",
#'     correction_for_sample_size = FALSE,
#'     letters = TRUE, add_points = TRUE,
#'     plot_with_tuckey = FALSE
#'   )
#'   if (requireNamespace("patchwork")) {
#'     patchwork::wrap_plots(p2, guides = "collect")
#'   }
#'   # Artificially modify data_fungi to force alpha-diversity effect
#'   data_fungi_modif <- clean_pq(subset_samples_pq(data_fungi, !is.na(data_fungi@sam_data$Height)))
#'   data_fungi_modif@otu_table[data_fungi_modif@sam_data$Height == "High", ] <-
#'     data_fungi_modif@otu_table[data_fungi_modif@sam_data$Height == "High", ] +
#'     sample(c(rep(0, ntaxa(data_fungi_modif) / 2), rep(100, ntaxa(data_fungi_modif) / 2)))
#'   p3 <- hill_pq(data_fungi_modif, "Height",
#'     letters = TRUE, vioplot = TRUE,
#'     add_points = TRUE
#'   )
#' }
#' }
#' @seealso [psmelt_samples_pq()] and [ggbetween_pq()]
hill_pq <- function(physeq,
                    fact = NULL,
                    variable = NULL,
                    hill_scales = c(0, 1, 2),
                    color_fac = NA,
                    letters = FALSE,
                    add_points = FALSE,
                    add_info = TRUE,
                    kruskal_test = TRUE,
                    one_plot = FALSE,
                    plot_with_tuckey = TRUE,
                    correction_for_sample_size = TRUE,
                    na_remove = TRUE,
                    vioplot = FALSE) {
  if (!is.null(variable)) {
    if (!is.null(fact)) {
      stop(
        "You must set only one parameter of variable or fact. This 2
        parameters are strictly equivalent."
      )
    } else {
      variable_fac <- variable
    }
  } else {
    if (!is.null(fact)) {
      variable_fac <- fact
    } else {
      stop("You must set the parameter fact.")
    }
  }
  var <- sym(variable_fac)
  if (is.na(color_fac)) {
    color_fac <- sym(variable_fac)
  } else {
    color_fac <- sym(color_fac)
  }

  physeq <- taxa_as_rows(physeq)
  if (na_remove || letters) {
    physeq <- subset_samples_pq(physeq, !is.na(physeq@sam_data[[fact]]))
  }
  physeq@sam_data[[fact]] <- as.factor(physeq@sam_data[[fact]])


  otu_hill <-
    vegan::renyi(t(physeq)@otu_table,
      scales = hill_scales,
      hill = TRUE
    )
  colnames(otu_hill) <- paste0("Hill_", hill_scales)

  df_hill <- data.frame(otu_hill, physeq@sam_data)
  df_hill[, seq_along(hill_scales)] <-
    apply(df_hill[, seq_along(hill_scales)], 2, as.numeric)

  p_var <-
    hill_tuckey_pq(
      physeq,
      modality = variable_fac,
      hill_scales = hill_scales,
      correction_for_sample_size = correction_for_sample_size
    )
  p_list <- list()

  if (kruskal_test) {
    kt_res <- list()
    for (i in seq_along(hill_scales)) {
      kt_res[[i]] <- kruskal.test(df_hill[, paste0("Hill_", hill_scales[[i]])], df_hill[, fact])
    }
    if (sum(sapply(kt_res, function(x) {
      x$p.value > 0.05
    })) > 0) {
      message(paste0(sum(sapply(kt_res, function(x) {
        x$p.value > 0.05
      })), " out of ", length(kt_res), " Hill scales do not show any global trends with you factor ", fact, ". Tuckey HSD plot is not informative for those Hill scales. Letters are not printed for those Hill scales"))
    }
  }

  for (i in seq_along(hill_scales)) {
    if (vioplot) {
      p_list[[i]] <-
        ggplot(df_hill, aes(
          x = .data[[paste0("Hill_", hill_scales[[i]])]],
          y = !!var
        )) +
        geom_violin(aes(colour = as.factor(!!color_fac))) +
        labs(x = paste0("Hill_", hill_scales[[i]]))
    } else {
      p_list[[i]] <-
        ggplot(df_hill, aes(group = !!var, x = .data[[paste0("Hill_", hill_scales[[i]])]])) +
        geom_boxplot(outlier.size = 2, aes(colour = as.factor(!!color_fac), y = !!var)) +
        labs(x = paste0("Hill_", hill_scales[[i]]))
    }

    if (add_points) {
      p_list[[i]] <-
        p_list[[i]] + geom_jitter(aes(y = !!var, colour = as.factor(!!color_fac)), alpha = 0.5)
    }
    if (add_info) {
      subtitle_plot <- paste0(
        "Nb of samples: '",
        paste0(
          names(table(physeq@sam_data[[variable_fac]])),
          sep = "' : ",
          table(physeq@sam_data[[variable_fac]]),
          collapse = " - '"
        )
      )
      if (kruskal_test) {
        subtitle_plot <- paste0(
          subtitle_plot, "\n",
          paste0(
            " Hill ", hill_scales[[i]],
            " -- Kruskal-Wallis chi-squared =",
            round(kt_res[[i]]$statistic, 2),
            "; df = ", kt_res[[i]]$parameter,
            "; p.value =", format.pval(kt_res[[i]]$p.value, 2)
          )
        )
      }
      p_list[[i]] <- p_list[[i]] + labs(subtitle = subtitle_plot)
    }

    if (letters) {
      data_h <-
        p_var$data[grep(paste0("Hill_", hill_scales[[i]]), p_var$data[, 5]), ]
      data_h_pval <- data_h$`p adj`
      names(data_h_pval) <- data_h$modality
      Letters <-
        multcompView::multcompLetters(data_h_pval, reversed = TRUE)$Letters

      dt <- data.frame(variab = names(Letters), Letters = Letters)
      names(dt) <- c(var, "Letters")
      data_letters <- p_list[[i]]$data %>%
        group_by(!!var) %>%
        summarize(pos_letters = max(.data[[paste0("Hill_", hill_scales[[i]])]]) + 1) %>%
        inner_join(dt, by = join_by(!!fact))

      if (!kruskal_test | kt_res[[i]]$p.value < 0.05) {
        p_list[[i]] <- p_list[[i]] +
          geom_label(
            data = data_letters,
            aes(
              x = pos_letters,
              label = Letters,
            ),
            y = unique(ggplot_build(p_list[[i]])$data[[1]]$y),
            size = 4,
            stat = "unique",
            parse = TRUE
          )
      }
    }
  }

  res <- p_list
  if (plot_with_tuckey) {
    res[["tuckey"]] <- p_var
  }

  if (one_plot) {
    requireNamespace("patchwork", quietly = TRUE)
    if (letters || !plot_with_tuckey) {
      res[["tuckey"]] <- NULL
    }
    res <- patchwork::wrap_plots(res)
  }
  return(res)
}
################################################################################

################################################################################
#' Box/Violin plots for between-subjects comparisons of Hill Number
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Note that contrary to [hill_pq()], this function does not take into
#' account for difference in the number of sequences per samples/modalities.
#' You may use rarefy_by_sample = TRUE if the mean number of sequences per
#' samples differs among modalities.
#'
#' Basically a wrapper of function [ggstatsplot::ggbetweenstats()] for
#' object of class phyloseq
#' @inheritParams clean_pq
#' @param fact (required): The variable to test. Must be present in
#'   the `sam_data` slot of the physeq object.
#' @param one_plot (logical, default FALSE) If TRUE, return a unique
#'   plot with the three plot inside using the patchwork package.
#' @param rarefy_by_sample (logical, default FALSE) If TRUE, rarefy
#'   samples using [phyloseq::rarefy_even_depth()] function
#' @param ... Other arguments passed on to [ggstatsplot::ggbetweenstats()] function.

#' @return Either an unique ggplot2 object (if one_plot is TRUE) or
#'  a list of 3 ggplot2 plot:
#' - plot_Hill_0 : the ggbetweenstats of Hill number 0 (= species richness)
#'     against the variable fact
#' - plot_Hill_1 : the ggbetweenstats of Hill number 1 (= Shannon index)
#'      against the variable fact
#' - plot_Hill_2 : the ggbetweenstats of Hill number 2 (= Simpson index)
#'     against the variable fact
#'
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("ggstatsplot")) {
#'   p <- ggbetween_pq(data_fungi, fact = "Time", p.adjust.method = "BH")
#'   p[[1]]
#'   ggbetween_pq(data_fungi, fact = "Height", one_plot = TRUE)
#'   ggbetween_pq(data_fungi, fact = "Height", one_plot = TRUE, rarefy_by_sample = TRUE)
#' }
#' }
#' @author Adrien Taudière
#' @details This function is mainly a wrapper of the work of others.
#'   Please make a reference to `ggstatsplot::ggbetweenstats()` if you
#'   use this function.

ggbetween_pq <-
  function(physeq,
           fact,
           one_plot = FALSE,
           rarefy_by_sample = FALSE,
           ...) {
    verify_pq(physeq)
    physeq <- taxa_as_columns(physeq)

    if (rarefy_by_sample) {
      physeq <- clean_pq(rarefy_even_depth(physeq))
    }

    if (are_modality_even_depth(physeq, fact)$p.value < 0.05) {
      warning(
        paste0(
          "The mean number of sequences per samples vary across modalities of the variable '",
          fact,
          "' You should use rarefy_by_sample = TRUE or try hill_pq() with correction_for_sample_size = TRUE"
        )
      )
    }

    df <- cbind(
      "nb_taxa" = sample_sums(physeq@otu_table),
      physeq@sam_data,
      "hill_0" = vegan::renyi(physeq@otu_table, scales = 0, hill = TRUE),
      "hill_1" = vegan::renyi(physeq@otu_table, scales = 1, hill = TRUE),
      "hill_2" = vegan::renyi(physeq@otu_table, scales = 2, hill = TRUE)
    )
    fact <- sym(fact)
    p0 <- ggstatsplot::ggbetweenstats(df, !!fact, hill_0, ...)
    p1 <- ggstatsplot::ggbetweenstats(df, !!fact, hill_1, ...)
    p2 <- ggstatsplot::ggbetweenstats(df, !!fact, hill_2, ...)

    res <- list(
      "plot_Hill_0" = p0,
      "plot_Hill_1" = p1,
      "plot_Hill_2" = p2
    )

    if (one_plot) {
      requireNamespace("patchwork", quietly = TRUE)
      res <- res[[1]] + res[[2]] + res[[3]]
    }
    return(res)
  }











################################################################################
#' Summarize a \code{\link{phyloseq-class}} object using a plot.
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Graphical representation of a phyloseq object.
#' 
#' @inheritParams clean_pq
#' @param add_info Does the bottom down corner contain
#'   extra informations?
#' @param min_seq_samples (int): Used only when add_info is set
#'   to true to print the number of samples with less sequences than
#'   this number.
#' @param clean_pq (logical): Does the phyloseq
#'   object is cleaned using the [clean_pq()] function?
#' @examples
#'
#' summary_plot_pq(data_fungi)
#' summary_plot_pq(data_fungi, add_info = FALSE) + scale_fill_viridis_d()
#' @return A ggplot2 object
#' @export
summary_plot_pq <- function(physeq,
                            add_info = TRUE,
                            min_seq_samples = 500,
                            clean_pq = TRUE) {
  if (clean_pq) {
    physeq <- clean_pq(physeq)
  }
  if (physeq@otu_table@taxa_are_rows) {
    otu_tab <- t(physeq@otu_table)
  } else {
    otu_tab <- physeq@otu_table
  }
  d <- data.frame(
    y1 = c(1, 1, 5.3, 1),
    y2 = c(5, 5, 7.5, 5),
    x1 = c(1, 3.15, 1, 4.3),
    x2 = c(3, 4.15, 3, 5.5),
    names = c("@otu_table", "@tax_table", "@sam_data", "@refseq"),
    nb_values = c(
      paste(
        format(ncol(otu_tab), big.mark = " "),
        "taxa\n",
        format(sum(otu_tab), big.mark = " "),
        "sequences\n",
        format(sum(otu_tab > 0), big.mark = " "),
        "occurrences"
      ),
      paste(ncol(physeq@tax_table), "taxonomic \n levels"),
      paste(
        ncol(physeq@sam_data),
        "variables\n",
        nsamples(physeq),
        "samples"
      ),
      paste(
        "Sequences length:\n",
        ifelse(
          is.null(physeq@refseq),
          "No refseq slot",
          paste(
            round(mean(
              Biostrings::width(physeq@refseq)
            ), 2),
            "+/-",
            round(stats::sd(
              Biostrings::width(physeq@refseq)
            ), 2)
          )
        )
      )
    )
  )

  p <- ggplot() +
    scale_x_continuous(name = "x") +
    scale_y_reverse(name = "y") +
    theme_void() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#aa4c26", "#003f5f", "khaki4", "#c8a734")) +
    geom_rect(
      data = d,
      mapping = aes(
        xmin = x1,
        xmax = x2,
        ymin = y1,
        ymax = y2,
        fill = names
      ),
      color = "black",
      alpha = 0.5
    ) +
    geom_text(
      data = d,
      aes(
        x = x1 + (x2 - x1) / 2,
        y = y1 + (y2 - y1) / 1.7,
        label = nb_values
      ),
      size = 4.5
    ) +
    geom_text(
      data = d,
      aes(
        x = x1 + (x2 - x1) / 2,
        y = y1 + (y2 - y1) / 5,
        label = names
      ),
      size = 6
    ) +
    annotate(
      geom = "text",
      x = 0.65,
      y = 3,
      label = "Taxa",
      size = 6,
      color = "#aa4c26",
      angle = 90,
      fontface = 2
    ) +
    annotate(
      geom = "text",
      x = 0.85,
      y = 3,
      label = "(OTUs, ASVs, ...)",
      size = 5,
      color = "#aa4c26",
      angle = 90
    ) +
    annotate(
      geom = "text",
      x = 2,
      y = 0.65,
      label = "Samples",
      size = 6,
      fontface = 2,
      color = "khaki4"
    )

  if (add_info) {
    supplementary_info <-
      data.frame(
        y1 = 5.3,
        y2 = 7.5,
        x1 = 3.15,
        nb_values =
          paste0(
            "Min nb seq per sample (",
            stringr::str_trunc(names(sort(
              sample_sums(otu_tab)
            ))[1], 15, "right"),
            min(sample_sums(otu_tab)),
            "\n",
            "Nb samples with less than ",
            min_seq_samples,
            " seq : ",
            sum(sample_sums(otu_tab) < min_seq_samples),
            "\n",
            "Min nb seq per taxa: ",
            min(taxa_sums(otu_tab)),
            " (",
            sum(taxa_sums(otu_tab) == min(taxa_sums(otu_tab))),
            " Taxons)",
            "\n",
            "Min seq length: ",
            ifelse(
              is.null(physeq@refseq),
              "No refseq slot",
              min(Biostrings::width(physeq@refseq))
            ),
            "\n",
            "Max nb seq 1 taxa in 1 sample: ",
            max(otu_tab),
            "\n",
            "Max nb of sample for one taxon (",
            names(sort(taxa_sums(otu_tab > 0), decreasing = TRUE))[1],
            "): ",
            max(taxa_sums(otu_tab > 0)),
            "\n",
            "Nb of taxa present in 1 sample only: ",
            sum(taxa_sums(otu_tab > 0) == 1)
          )
      )

    p <- p +
      geom_text(
        data = supplementary_info,
        aes(
          x = x1,
          y = y1 + (y2 - y1) / 2.1,
          label = nb_values
        ),
        size = 3.5,
        hjust = 0
      )
  }

  return(p)
}
################################################################################

################################################################################
#' rotl wrapper for phyloseq data
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Make a phylogenetic tree using the ASV names of a physeq object and the
#'   Open Tree of Life tree.
#'
#' @inheritParams clean_pq
#' @param species_colnames (default: "Genus_species"): the name of the column
#'   where the species binominal name is stored in `@tax_table` slot.
#'   Can also be a vector of two columns names e.g. `c("Genus", "Species")`
#' @param context_name : can bue used to select only a part of the Open Tree
#'   of Life. See `?rotl::tnrs_contexts()` for available values
#' @return A plot
#' @export
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `rotl` package if you
#'   use this function.
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' \donttest{
#' if (requireNamespace("rotl")) {
#'   tr <- rotl_pq(data_fungi_mini, species_colnames = "Genus_species")
#'   plot(tr)
#'
#'   tr_Asco <- rotl_pq(data_fungi, species_colnames = "Genus_species", context_name = "Ascomycetes")
#'   plot(tr_Asco)
#' }
#' }
rotl_pq <- function(physeq,
                    species_colnames = "Genus_species",
                    context_name = "All life") {
  if (length(species_colnames) == 2) {
    physeq@tax_table <- tax_table(cbind(
      physeq@tax_table,
      "Genus_species" = paste(physeq@tax_table[, species_colnames[1]],
        physeq@tax_table[, species_colnames[2]],
        sep = "_"
      )
    ))
    species_colnames <- "Genus_species"
  }
  taxa_names_rotl <- physeq@tax_table[, species_colnames]
  taxa_names_rotl <- taxa_names_rotl[!grepl("NA", taxa_names_rotl)]
  taxa_names_rotl <- c(unclass(gsub("_", " ", taxa_names_rotl)))

  resolved_names <- rotl::tnrs_match_names(taxa_names_rotl)
  resolved_names <- resolved_names[resolved_names$flags == "", ]
  clean_taxa_names_rotl <-
    taxa_names_rotl[taxa_names_rotl %in% resolved_names$unique_name]

  resolved_names2 <-
    rotl::tnrs_match_names(clean_taxa_names_rotl, context_name = context_name)

  tr <-
    rotl::tol_induced_subtree(ott_ids = rotl::ott_id(resolved_names2))
  return(tr)
}
################################################################################

################################################################################
#' Heat tree from `metacoder` package using `tax_table` slot
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Note that the number of ASV is store under the name `n_obs`
#' and the number of sequences under the name `nb_sequences`
#'
#' @inheritParams clean_pq
#' @param taxonomic_level (default: NULL): a vector of selected
#' taxonomic level using
#'   their column numbers (e.g. taxonomic_level = 1:7)
#' @param ... Arguments passed on to \code{\link[metacoder]{heat_tree}}
#'
#' @return A plot
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' \donttest{
#' if (requireNamespace("metacoder")) {
#'   library("metacoder")
#'   data("GlobalPatterns", package = "phyloseq")
#'
#'   GPsubset <- subset_taxa(
#'     GlobalPatterns,
#'     GlobalPatterns@tax_table[, 1] == "Bacteria"
#'   )
#'
#'   GPsubset <- subset_taxa(
#'     GPsubset,
#'     rowSums(GPsubset@otu_table) > 5000
#'   )
#'
#'   GPsubset <- subset_taxa(
#'     GPsubset,
#'     rowSums(is.na(GPsubset@tax_table)) == 0
#'   )
#'
#'   heat_tree_pq(GPsubset,
#'     node_size = n_obs,
#'     node_color = n_obs,
#'     node_label = taxon_names,
#'     tree_label = taxon_names,
#'     node_size_trans = "log10 area"
#'   )
#'
#'   heat_tree_pq(GPsubset,
#'     node_size = nb_sequences,
#'     node_color = n_obs,
#'     node_label = taxon_names,
#'     tree_label = taxon_names,
#'     node_size_trans = "log10 area"
#'   )
#' }
#' }
heat_tree_pq <- function(physeq, taxonomic_level = NULL, ...) {
  requireNamespace("metacoder", quietly = TRUE)
  if (!is.null(taxonomic_level)) {
    physeq@tax_table <- physeq@tax_table[, taxonomic_level]
  }

  data_metacoder <- metacoder::parse_phyloseq(physeq)
  data_metacoder$data$taxon_counts <-
    metacoder::calc_taxon_abund(data_metacoder, data = "otu_table")
  data_metacoder$data$taxon_counts$nb_sequences <-
    rowSums(data_metacoder$data$taxon_counts[, -1])

  p <- heat_tree(data_metacoder, ...)

  return(p)
}
################################################################################

################################################################################
#' Visualization of two samples for comparison
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Graphical representation of distribution of taxa across two samples.
#' 
#' @inheritParams clean_pq
#' @param fact (default: NULL) Name of the factor in `physeq@sam_data`.
#'   If left to NULL use the `left_name` and `right_name` parameter as modality.
#' @param merge_sample_by (default: NULL) if not `NULL` samples of
#'   physeq are merged using the vector set by `merge_sample_by`. This
#'   merging used the [merge_samples2()]. In the case of
#'   [biplot_pq()] this must be a factor with two levels only.
#' @param rarefy_after_merging Rarefy each sample after merging by the
#'   modalities merge_sample_by
#' @param inverse_side Inverse the side (put the right modality in the left side).
#' @param left_name Name fo the left sample.
#' @param left_name_col Color for the left name
#' @param left_fill Fill fo the left sample.
#' @param left_col Color fo the left sample.
#' @param right_name Name fo the right sample.
#' @param right_name_col Color for the right name
#' @param right_fill Fill fo the right sample.
#' @param right_col Color fo the right sample.
#' @param log10trans (logical) Does abundancy is log10 transformed ?
#' @param nudge_y A parameter to control the y position of abundancy values.
#'   If a vector of two values are set. The first value is for the left side.
#'   and the second value for the right one. If one value is set,
#'   this value is used for both side.
#' @param geom_label (default: FALSE, logical) if TRUE use the [ggplot2::geom_label()] function
#'   instead of [ggplot2::geom_text()] to indicate the numbers of sequences.
#' @param text_size size for the number of sequences
#' @param size_names size for the names of the 2 samples
#' @param y_names y position for the names of the 2 samples. If NA (default),
#'   computed using the maximum abundances values.
#' @param ylim_modif vector of two values. Modificator (by a multiplication)
#'   of ylim. If one value is set, this value is used for both limits.
#' @param nb_samples_info (default: TRUE, logical) if TRUE and merge_sample_by is set,
#'   add the number of samples merged for both levels.
#' @param plotly_version If TRUE, use [plotly::ggplotly()] to return
#'   a interactive ggplot.
#' @param ... Other arguments for the ggplot function
#' @return A plot
#'
#' @examples
#' data_fungi_2Height <- subset_samples(data_fungi_mini, Height %in% c("Low", "High"))
#' biplot_pq(data_fungi_2Height, "Height", merge_sample_by = "Height")
#' @export
#' @author Adrien Taudière
#'
biplot_pq <- function(physeq,
                      fact = NULL,
                      merge_sample_by = NULL,
                      rarefy_after_merging = FALSE,
                      inverse_side = FALSE,
                      left_name = NULL,
                      left_name_col = "#4B3E1E",
                      left_fill = "#4B3E1E",
                      left_col = "#f3f2d9",
                      right_name = NULL,
                      right_name_col = "#1d2949",
                      right_fill = "#1d2949",
                      right_col = "#1d2949",
                      log10trans = TRUE,
                      nudge_y = c(0.3, 0.3),
                      geom_label = FALSE,
                      text_size = 3,
                      size_names = 5,
                      y_names = NA,
                      ylim_modif = c(1, 1),
                      nb_samples_info = TRUE,
                      plotly_version = FALSE,
                      ...) {
  if (!is.null(merge_sample_by)) {
    if (nb_samples_info) {
      modality_1_nb <- table(physeq@sam_data[, merge_sample_by])[1]
      modality_2_nb <- table(physeq@sam_data[, merge_sample_by])[2]
    }
    physeq <- merge_samples2(physeq, merge_sample_by)
    physeq <- clean_pq(physeq)
  }

  if (nsamples(physeq) != 2) {
    stop(
      "biplot_pq needs only two samples in the
    physeq object or a valid merge_sample_by parameter"
    )
  }

  if (rarefy_after_merging) {
    physeq <- clean_pq(rarefy_even_depth(physeq))
  }

  if (sample_sums(physeq)[1] / sample_sums(physeq)[2] > 2 ||
    sample_sums(physeq)[2] / sample_sums(physeq)[1] > 2) {
    message(
      paste0(
        "The two modalities differ greatly (more than x2) in their number of sequences (",
        sample_sums(physeq)[1],
        " vs ",
        sample_sums(physeq)[2],
        "). You may be interested by the parameter rarefy_after_merging"
      )
    )
  }

  if (is.null(fact)) {
    if (is.null(left_name)) {
      left_name <- "A"
    }
    if (is.null(right_name)) {
      right_name <- "B"
    }
    modality <- factor(c(left_name, right_name))
  } else {
    modality <-
      as.factor(eval(parse(
        text = paste("physeq@sam_data$", fact, sep = "")
      )))
  }

  if (inverse_side) {
    modality <- factor(modality, rev(levels(as.factor(modality))))
  }

  if (is.null(left_name)) {
    left_name <- levels(modality)[1]
  }
  if (is.null(right_name)) {
    right_name <- levels(modality)[2]
  }

  if (!is.null(merge_sample_by) && nb_samples_info) {
    left_name <- paste0(left_name, " (", modality_1_nb, " samples)")
    right_name <-
      paste0(right_name, " (", modality_2_nb, " samples)")
  }

  physeq@sam_data$modality <- modality

  mdf <- phyloseq::psmelt(physeq)
  mdf <- mdf[mdf$Abundance > 0, ]
  # mdf <- dplyr::rename(mdf, Abundance = Abundance)

  if (length(ylim_modif) == 1) {
    ylim_modif <- c(ylim_modif, ylim_modif)
  }

  if (length(y_names) == 1) {
    y_names <- c(y_names, y_names)
  }

  if (length(nudge_y) == 1) {
    nudge_y <- c(nudge_y, nudge_y)
  }
  if (log10trans) {
    mdf$Ab <- log10(mdf$Abundance + 1)
  } else {
    mdf$Ab <- mdf$Abundance
    nudge_y <- mean(mdf$Abundance) * nudge_y
  }

  mdf$Ab[mdf$modality == levels(modality)[1]] <-
    -mdf$Ab[mdf$modality == levels(modality)[1]]
  mdf$Proportion <- paste0(round(100 * mdf$Abundance /
    sum(mdf$Abundance[mdf$modality == levels(modality)[2]]), 2), "%")
  mdf$Proportion[mdf$modality == levels(modality)[1]] <-
    paste0(round(100 * mdf$Abundance[mdf$modality == levels(modality)[1]] /
      sum(mdf$Abundance[mdf$modality == levels(modality)[1]]), 2), "%")

  p <- mdf %>%
    ggplot(
      aes(
        x = stats::reorder(OTU, Abundance),
        y = Ab,
        fill = modality,
        names = OTU,
        Ab = Abundance,
        Proportion = Proportion,
        Family = Family,
        Genus = Genus,
        Species = Species
      ),
      ...
    ) +
    geom_bar(stat = "identity", width = .6) +
    annotate(
      "rect",
      xmin = "Samples",
      xmax = "Samples",
      ymin = -max(mdf$Ab),
      ymax = max(mdf$Ab)
    ) +
    annotate(
      geom = "text",
      label = right_name,
      x = "Samples",
      y = ifelse(is.na(y_names), max(mdf$Ab) / 2, y_names[2]),
      hjust = 0.5,
      vjust = 0.5,
      size = size_names,
      fontface = "bold",
      col = right_name_col
    ) +
    annotate(
      geom = "text",
      label = left_name,
      x = "Samples",
      y = ifelse(is.na(y_names), (min(mdf$Ab) / 2), -y_names[1]),
      hjust = 0.5,
      vjust = 0.5,
      size = size_names,
      fontface = "bold",
      col = left_name_col
    ) +
    geom_hline(aes(yintercept = 0)) +
    scale_x_discrete(limits = c(names(sort(
      tapply(mdf$Abundance, mdf$OTU, sum)
    )), "Samples")) +
    ylim(min(mdf$Ab), max(mdf$Ab) * 1.1)


  if (geom_label) {
    p <- p +
      geom_label(
        aes(
          label = Abundance,
          color = modality,
          fill = modality,
          alpha = 0.5,
          y = ifelse(Ab > 0, Ab + nudge_y[2], Ab + nudge_y[1])
        ),
        size = text_size
      )
  } else {
    p <- p +
      geom_text(
        aes(
          label = Abundance,
          color = modality,
          y = ifelse(Ab > 0, Ab + nudge_y[2], Ab + nudge_y[1])
        ),
        size = text_size
      )
  }

  p <- p + coord_flip() +
    theme_minimal() +
    theme(plot.title = element_text(hjust = .5), axis.ticks = element_blank()) +
    scale_fill_manual(values = c(left_fill, right_fill)) +
    scale_color_manual(values = c(left_col, right_col), guide = "none") +
    ylim(
      c(
        layer_scales(p)$y$get_limits()[1] * ylim_modif[1],
        layer_scales(p)$y$get_limits()[2] * ylim_modif[2]
      )
    )

  if (plotly_version) {
    p <- plotly::ggplotly(
      p,
      tooltip = c("OTU", "Ab", "Proportion", "Family", "Genus", "Species"),
      height = 1200,
      width = 800
    ) %>%
      plotly::layout(
        xaxis = list(autorange = TRUE),
        yaxis = list(autorange = TRUE)
      ) %>%
      plotly::config(locale = "fr") %>%
      plotly::hide_legend()
  }
  return(p)
}
################################################################################


################################################################################
#' Visualization of a collection of couples of samples for comparison
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' This allow to plot all the possible [biplot_pq()] combination
#' using one factor.
#'
#' @inheritParams clean_pq
#' @param split_by (required if pairs is NULL) the name of the factor to make all combination
#'   of couples of values
#' @param pairs (required if pairs is NULL) the name of the factor in physeq@sam_data` slot
#'   to make plot by pairs of samples. Each level must be present only two times.
#'   Note that if you set pairs, you also must set fact arguments to pass on to [biplot_pq()].
#' @param na_remove (logical, default TRUE) if TRUE remove all the samples
#'   with NA in the `split_by` variable of the `physeq@sam_data` slot
#' @param ... Other parameters passed on to [biplot_pq()]
#'
#' @return a list of ggplot object
#' @export
#'
#' @examples
#' \donttest{
#' data_fungi_abun <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)
#' p <- multi_biplot_pq(data_fungi_abun, "Height")
#' lapply(p, print)
#' }
#' @author Adrien Taudière
multi_biplot_pq <- function(physeq,
                            split_by = NULL,
                            pairs = NULL,
                            na_remove = TRUE,
                            ...) {
  if (is.null(pairs) && is.null(split_by)) {
    stop("You must set one of split_by or pairs.")
  } else if (!is.null(pairs) && !is.null(split_by)) {
    stop("You must set either split_by or pairs, not both.")
  } else if (!is.null(split_by) &&
    is.null(physeq@sam_data[[split_by]])) {
    stop("split_by must be set and must be a variable in physeq@sam_data")
  } else if (!is.null(pairs) && is.null(physeq@sam_data[[pairs]])) {
    stop("pairs must be set and must be a variable in physeq@sam_data")
  }

  if (na_remove && !is.null(split_by)) {
    new_physeq <-
      subset_samples_pq(physeq, !is.na(physeq@sam_data[[split_by]]))
    if (nsamples(physeq) - nsamples(new_physeq) > 0) {
      message(
        paste0(
          nsamples(physeq) - nsamples(new_physeq),
          " were discarded due to NA in variables present in formula."
        )
      )
    }
    physeq <- new_physeq
  }

  if (!is.null(pairs)) {
    p <- list()
    for (c in levels(as.factor(physeq@sam_data[[pairs]]))) {
      new_physeq <-
        subset_samples_pq(physeq, physeq@sam_data[[pairs]] %in% c)
      p[[c]] <- biplot_pq(new_physeq, ...) + ggtitle(c)
    }
  } else {
    names_split_by <- names(table(physeq@sam_data[[split_by]]))
    couples <- combn(names_split_by, 2)

    p <- list()
    for (c in seq_along(ncol(couples))) {
      names_p <- paste0(couples[1, c], " - ", couples[2, c])
      new_physeq <-
        subset_samples_pq(physeq, physeq@sam_data[[split_by]] %in%
          c(couples[1, c], couples[2, c]))
      p[[names_p]] <- biplot_pq(new_physeq,
        fact = split_by,
        merge_sample_by = split_by,
        ...
      )
    }
  }
  return(p)
}
################################################################################


################################################################################
#' Plot taxonomic distribution in function of a factor with stacked bar in %
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   An alternative to `phyloseq::plot_bar()` function.
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}.
#' @param merge_sample_by a vector to determine
#'   which samples to merge using the
#'   [merge_samples2()] function.
#'   Need to be in \code{physeq@sam_data}
#' @param type If "nb_seq" (default), the number of sequences is
#'   used in plot. If "nb_taxa", the number of ASV is plotted. If both,
#'   return a list of two plots, one for nbSeq and one for ASV.
#' @param taxa_fill (default: 'Order'): Name of the taxonomic rank of interest
#' @param print_values (logical, default TRUE): Do we print some values on plot?
#' @param color_border color for the border
#' @param linewidth The line width of geom_bar
#' @param prop_print_value minimal proportion to print value (default 0.01)
#' @param nb_print_value number of higher values to print
#'    (replace prop_print_value if both are set).
#' @param add_info (logical, default TRUE) Do we add title and subtitle with
#'   information about the total number of sequences and the number of samples
#'   per modality.
#' @param na_remove (logical, default TRUE) if TRUE remove all the samples
#'   with NA in the `split_by` variable of the `physeq@sam_data` slot
#' @param clean_pq (logical)
#'   If set to TRUE, empty samples are discarded after subsetting ASV
#' @return A ggplot2 graphic
#' @export
#' @author Adrien Taudière
#' @seealso [tax_bar_pq()] and [multitax_bar_pq()]
#' @examples
#' data(data_fungi_sp_known)
#' plot_tax_pq(data_fungi_sp_known,
#'   "Time",
#'   merge_sample_by = "Time",
#'   taxa_fill = "Class"
#' )
#' \donttest{
#' plot_tax_pq(data_fungi_sp_known,
#'   "Height",
#'   merge_sample_by = "Height",
#'   taxa_fill = "Class",
#'   na_remove = TRUE,
#'   color_border = rgb(0, 0, 0, 0)
#' )
#'
#' plot_tax_pq(data_fungi_sp_known,
#'   "Height",
#'   merge_sample_by = "Height",
#'   taxa_fill = "Class",
#'   na_remove = FALSE,
#'   clean_pq = FALSE
#' )
#' }
plot_tax_pq <-
  function(physeq,
           fact = NULL,
           merge_sample_by = NULL,
           type = "nb_seq",
           taxa_fill = "Order",
           print_values = TRUE,
           color_border = "lightgrey",
           linewidth = 0.1,
           prop_print_value = 0.01,
           nb_print_value = NULL,
           add_info = TRUE,
           na_remove = TRUE,
           clean_pq = TRUE) {
    if (na_remove) {
      new_physeq <-
        subset_samples_pq(physeq, !is.na(physeq@sam_data[[fact]]))
      if (nsamples(physeq) - nsamples(new_physeq) > 0) {
        message(
          paste0(
            nsamples(physeq) - nsamples(new_physeq),
            " were discarded due to NA in variables present in formula."
          )
        )
      }
      physeq <- new_physeq
    }

    if (clean_pq) {
      physeq <- clean_pq(physeq)
    }


    physeq_old <- physeq

    if (!is.null(merge_sample_by)) {
      physeq <- merge_samples2(physeq, merge_sample_by)
    }

    if (!is.null(nb_print_value)) {
      prop_print_value <-
        taxa_sums(physeq)[nb_print_value] / sum(physeq@otu_table)
    }

    if (type %in% c("nb_seq", "both")) {
      mdf <- psmelt(physeq)
      mdf <- mdf %>% mutate(percent = Abundance / sum(Abundance))

      p_seq <-
        ggplot(mdf, aes(x = .data[[fact]], y = .data[["Abundance"]], fill = .data[[taxa_fill]])) +
        geom_bar(
          aes(fill = .data[[taxa_fill]]),
          stat = "identity",
          position = "fill",
          color = color_border,
          linewidth = linewidth
        ) +
        scale_y_continuous(labels = scales::percent) +
        ylab("Pseudo-abundance (nb of sequences)")

      if (print_values) {
        p_seq <- p_seq + geom_text(
          aes(label = round(Abundance)),
          color = ifelse(
            p_seq$data$percent > prop_print_value,
            "black",
            rgb(1, 1, 1, 0)
          ),
          position = position_fill(vjust = 0.5)
        )
      }
    }
    if (type %in% c("nb_taxa", "both")) {
      mdf <-
        psmelt(as_binary_otu_table(physeq))
      mdf <- mdf %>% mutate(percent = Abundance / sum(Abundance))

      p_taxa <-
        ggplot(mdf, aes(x = .data[[fact]], y = .data[["Abundance"]], fill = .data[[taxa_fill]])) +
        geom_bar(
          aes(fill = .data[[taxa_fill]]),
          stat = "identity",
          position = "fill",
          color = color_border,
          linewidth = linewidth
        ) +
        scale_y_continuous(labels = scales::percent) +
        ylab("Nb_ASV")
    }

    if (add_info) {
      if (type %in% c("nb_seq", "both")) {
        p_seq <- p_seq +
          labs(
            title = paste("Total nb of sequences: ", sum(physeq_old@otu_table)),
            subtitle = paste0(
              "Nb of samples: '",
              paste0(
                names(table(physeq_old@sam_data[[fact]])),
                sep = "' : ",
                table(physeq_old@sam_data[[fact]]),
                collapse = " - '"
              )
            )
          )
      }
      if (type %in% c("nb_taxa", "both")) {
        p_taxa <- p_taxa +
          labs(
            title = paste("Total nb of sequences: ", sum(physeq_old@otu_table)),
            subtitle = paste0(
              "Nb of samples: '",
              paste0(
                names(table(physeq_old@sam_data[[fact]])),
                sep = "' : ",
                table(physeq_old@sam_data[[fact]]),
                collapse = " - '"
              )
            )
          )
      }
    }

    if (type == "nb_seq") {
      return(p_seq)
    } else if (type == "nb_taxa") {
      return(p_taxa)
    } else if (type == "both") {
      return(list(p_seq, p_taxa))
    }
  }
################################################################################


################################################################################
#' Plot taxonomic distribution across 3 taxonomic levels and optionally
#' one sample factor
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Note that lvl3 need to be nested in lvl2 which need to be nested
#' in lvl1
#'
#' @inheritParams clean_pq
#' @param lvl1 (required) Name of the first (higher) taxonomic rank of interest
#' @param lvl2 (required) Name of the second (middle) taxonomic rank of interest
#' @param lvl3 (required) Name of the first (lower) taxonomic rank of interest
#' @param fact Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}. If not set, the taxonomic
#'   distribution is plot for all samples together.
#' @param nb_seq (logical; default TRUE) If set to FALSE, only the number of ASV
#'   is count. Concretely, physeq otu_table is transformed in a binary
#'   otu_table (each value different from zero is set to one)
#' @param log10trans (logical, default TRUE) If TRUE,
#'   the number of sequences (or ASV if nb_seq = FALSE) is log10
#'   transformed.
#' @return A ggplot2 graphic
#' @export
#'
#' @author Adrien Taudière
#' @examples
#' \donttest{
#' if (requireNamespace("ggh4x")) {
#'   multitax_bar_pq(data_fungi_sp_known, "Phylum", "Class", "Order", "Time")
#'   multitax_bar_pq(data_fungi_sp_known, "Phylum", "Class", "Order")
#'   multitax_bar_pq(data_fungi_sp_known, "Phylum", "Class", "Order",
#'     nb_seq = FALSE, log10trans = FALSE
#'   )
#' }
#' }
multitax_bar_pq <- function(physeq,
                            lvl1,
                            lvl2,
                            lvl3,
                            fact = NULL,
                            nb_seq = TRUE,
                            log10trans = TRUE) {
  psm_1 <- psmelt(physeq) %>%
    filter(Abundance > 0) %>%
    filter(!is.na(.data[[lvl1]])) %>%
    filter(!is.na(.data[[lvl3]])) %>%
    filter(!is.na(.data[[lvl2]]))

  if (is.null(fact)) {
    psm_2 <- psm_1 %>%
      group_by(OTU) %>%
      summarise(Abundance = sum(Abundance))

    psm <- inner_join(psm_2,
      psm_1[, c("OTU", lvl1, lvl2, lvl3)],
      by = join_by("OTU" == "OTU"),
      multiple =
        "first"
    )

    if (!nb_seq) {
      psm$Abundance <- 1
    }

    data_gg <- tibble(
      "Abundance" = tapply(psm$Abundance, psm[[lvl3]], sum),
      "LVL1" = tapply(psm[[lvl1]], psm[[lvl3]], unique),
      "LVL2" = tapply(psm[[lvl2]], psm[[lvl3]], unique),
      "LVL3" = tapply(psm[[lvl3]], psm[[lvl3]], unique)
    )

    if (log10trans) {
      data_gg$Abundance <- log10(data_gg$Abundance)
    }

    p <- ggplot(data_gg, aes(
      x = Abundance,
      fill = LVL1,
      y = LVL3
    )) +
      geom_bar(stat = "identity") +
      ggh4x::facet_nested(LVL1 + LVL2 ~ ., scales = "free", space = "free") +
      theme(strip.text.y.right = element_text(angle = 0)) +
      theme(legend.position = "none")
  } else {
    psm_2 <- psm_1 %>%
      group_by(OTU, .data[[fact]]) %>%
      summarise(Abundance = sum(Abundance)) %>%
      filter(Abundance > 0)

    psm <- inner_join(psm_2,
      psm_1[, c("OTU", lvl1, lvl2, lvl3)],
      by = join_by("OTU" == "OTU"),
      multiple =
        "first"
    )

    if (!nb_seq) {
      psm$Abundance <- 1
    }

    data_gg <- tibble(
      "Abundance" = tapply(psm$Abundance, paste(psm[[fact]], psm[[lvl3]]), sum),
      "FACT" = tapply(psm[[fact]], paste(psm[[fact]], psm[[lvl3]]), unique),
      "LVL1" = tapply(psm[[lvl1]], paste(psm[[fact]], psm[[lvl3]]), unique),
      "LVL2" = tapply(psm[[lvl2]], paste(psm[[fact]], psm[[lvl3]]), unique),
      "LVL3" = tapply(psm[[lvl3]], paste(psm[[fact]], psm[[lvl3]]), unique)
    )

    if (log10trans) {
      data_gg$Abundance <- log10(data_gg$Abundance)
    }

    p <- ggplot(data_gg, aes(
      x = Abundance,
      fill = LVL1,
      y = LVL3
    )) +
      geom_bar(stat = "identity") +
      ggh4x::facet_nested(LVL1 + LVL2 ~ FACT, scales = "free", space = "free") +
      theme(strip.text.y.right = element_text(angle = 0)) +
      theme(legend.position = "none")
  }
  return(p)
}
################################################################################


################################################################################
#' Compute tSNE position of samples from a phyloseq object
#'
#' @inheritParams clean_pq
#' @param method A method to calculate distance using `vegan::vegdist()` function
#' @param dims (Int) Output dimensionality (default: 2)
#' @param theta (Numeric) Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.0 see details in the man page of `Rtsne::Rtsne`).
#' @param perplexity (Numeric) Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details in the man page of `Rtsne::Rtsne`)
#' @param ... Other arguments passed on to `Rtsne::Rtsne()`
#'
#' @return A list of element including the matrix Y containing the new representations for the objects.
#'   See ?Rtsne::Rtsne() for more information
#' @export
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' if (requireNamespace("Rtsne")) {
#'   res_tsne <- tsne_pq(data_fungi)
#' }
tsne_pq <-
  function(physeq,
           method = "bray",
           dims = 2,
           theta = 0.0,
           perplexity = 30,
           ...) {
    physeq <- taxa_as_rows(physeq)

    res_tsne <-
      Rtsne::Rtsne(
        vegan::vegdist(t(physeq@otu_table), method = method),
        dims = dims,
        theta = theta,
        perplexity = perplexity,
        is_distance = TRUE,
        ...
      )

    return(res_tsne)
  }
################################################################################


################################################################################
#' Plot a tsne low dimensional representation of a phyloseq object
#'
#' @description
#' 
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Partially inspired by `phylosmith::tsne_phyloseq()` function developed by Schuyler D. Smith.
#'
#' @inheritParams clean_pq
#' @param method A method to calculate distance using `vegan::vegdist()` function (default: "bray")
#' @param dims (Int) Output dimensionality (default: 2)
#' @param theta (Numeric) Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact TSNE (default: 0.0 see details in the man page of `Rtsne::Rtsne`).
#' @param perplexity (Numeric) Perplexity parameter (should not be bigger than 3 * perplexity < nrow(X) - 1, see details in the man page of `Rtsne::Rtsne`)
#' @param fact Name of the column in `physeq@sam_data` used to color points and compute ellipses.
#' @param ellipse_level The level used in stat_ellipse. Set to NULL to discard ellipse (default = 0.95)
#' @param plot_dims A vector of 2 values defining the rank of dimension to plot (default: c(1,2))
#' @param na_remove (logical, default TRUE) Does the samples with NA values in fact are removed? (default: true)
#' @param force_factor (logical, default TRUE) Force the fact column to be a factor.
#' @param ... Other arguments passed on to `Rtsne::Rtsne()`
#'
#' @return
#' A ggplot object
#'
#' @export
#' @author Adrien Taudière
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' data(data_fungi)
#' if (requireNamespace("Rtsne")) {
#'   plot_tsne_pq(data_fungi, fact = "Height", perplexity = 15)
#' }
#' \donttest{
#' if (requireNamespace("Rtsne")) {
#'   plot_tsne_pq(data_fungi, fact = "Time") + geom_label(aes(label = Sample_id, fill = Time))
#'   plot_tsne_pq(data_fungi, fact = "Time", na_remove = FALSE, force_factor = FALSE)
#' }
#' }
#'
plot_tsne_pq <- function(physeq,
                         method = "bray",
                         dims = 2,
                         theta = 0.0,
                         perplexity = 30,
                         fact = NA,
                         ellipse_level = 0.95,
                         plot_dims = c(1, 2),
                         na_remove = TRUE,
                         force_factor = TRUE,
                         ...) {
  if (!is.factor(physeq@sam_data[[fact]]) &&
    !is.na(fact) && force_factor) {
    physeq@sam_data[[fact]] <- as.factor(physeq@sam_data[[fact]])
  }

  if (na_remove && !is.na(fact)) {
    physeq <- subset_samples_pq(physeq, !is.na(physeq@sam_data[[fact]]))
  }

  tsne <- tsne_pq(
    physeq = physeq,
    method = method,
    dims = dims,
    theta = theta,
    perplexity = perplexity,
    ...
  )

  res_tSNE_A <- tsne$Y[, plot_dims[1]] / 100
  res_tSNE_B <- tsne$Y[, plot_dims[2]] / 100

  df <- tibble(
    res_tSNE_A,
    res_tSNE_B,
    as(physeq@sam_data, "data.frame")
  )


  g <-
    ggplot(data = df, aes(.data[["res_tSNE_A"]], .data[["res_tSNE_B"]], group = .data[[fact]])) +
    xlab(paste0("Dimension ", plot_dims[1], " of tSNE analysis")) +
    ylab(paste0("Dimension ", plot_dims[2], " of tSNE analysis"))

  g <- g + geom_point(
    aes(fill = .data[[fact]]),
    shape = 21,
    color = "black",
    size = 3,
    alpha = 1.0
  )

  if (!is.null(ellipse_level) && !is.na(fact)) {
    g <-
      g + stat_ellipse(aes(color = .data[[fact]]), level = ellipse_level)
  }

  return(g)
}
################################################################################


################################################################################
#' Scaling with ranked subsampling (SRS) curve of phyloseq object
#' 
#' @description 
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' A wraper of [SRS::SRScurve()] function.
#' @inheritParams clean_pq
#' @param clean_pq (logical): Does the phyloseq
#'   object is cleaned using the [clean_pq()] function?
#' @param ... Other arguments passed on to `SRS::SRScurve()`
#' @return A plot
#' @export
#'
#' @examples
#' if (requireNamespace("SRS")) {
#'   SRS_curve_pq(data_fungi_mini,
#'     max.sample.size = 200,
#'     rarefy.comparison = TRUE, rarefy.repeats = 3
#'   )
#'   SRS_curve_pq(data_fungi_mini, max.sample.size = 500, metric = "shannon")
#' }
SRS_curve_pq <- function(physeq, clean_pq = FALSE, ...) {
  if (clean_pq) {
    physeq <- clean_pq(physeq)
  }

  physeq <- taxa_as_rows(physeq)

  df <- data.frame(physeq@otu_table)

  SRS::SRScurve(df, ...)
}
################################################################################


################################################################################
#' iNterpolation and EXTrapolation of Hill numbers (with iNEXT)
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Note that this function is quite time-consuming due to high dimensionality in metabarcoding community matrix.
#' 
#' @inheritParams clean_pq
#' @param merge_sample_by (default: NULL) if not `NULL` samples of
#'   physeq are merged using the vector set by `merge_sample_by`. This
#'   merging used the [merge_samples2()]. In the case of
#'   [biplot_pq()] this must be a factor with two levels only.
#' @param ... Other arguments for the [iNEXT::iNEXT()] function
#' @return see [iNEXT::iNEXT()] documentation
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("iNEXT")) {
#'   data("GlobalPatterns", package = "phyloseq")
#'   GPsubset <- subset_taxa(
#'     GlobalPatterns,
#'     GlobalPatterns@tax_table[, 1] == "Bacteria"
#'   )
#'   GPsubset <- subset_taxa(
#'     GPsubset,
#'     rowSums(GPsubset@otu_table) > 20000
#'   )
#'   GPsubset <- subset_taxa(
#'     GPsubset,
#'     rowSums(is.na(GPsubset@tax_table)) == 0
#'   )
#'   GPsubset@sam_data$human <- GPsubset@sam_data$SampleType %in%
#'     c("Skin", "Feces", "Tong")
#'   res_iNEXT <- iNEXT_pq(
#'     GPsubset,
#'     merge_sample_by = "human",
#'     q = 1,
#'     datatype = "abundance",
#'     nboot = 2
#'   )
#'   iNEXT::ggiNEXT(res_iNEXT)
#'   iNEXT::ggiNEXT(res_iNEXT, type = 2)
#'   iNEXT::ggiNEXT(res_iNEXT, type = 3)
#' }
#' }
#' @author Adrien Taudière
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `iNEXT::iNEXT()` if you
#'   use this function.
#'
iNEXT_pq <- function(physeq,
                     merge_sample_by = NULL,
                     ...) {
  if (!is.null(merge_sample_by)) {
    physeq <- merge_samples2(physeq, merge_sample_by)
    physeq <- taxa_as_columns(physeq)
  }

  df <- data.frame(t(as.matrix(unclass(physeq@otu_table))))
  res_iNEXT <- iNEXT::iNEXT(df, ...)
  return(res_iNEXT)
}
################################################################################


################################################################################
#' Make upset plot for phyloseq object.
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Alternative to venn plot.
#'
#' @inheritParams clean_pq
#' @param fact (required): Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}.
#' @param min_nb_seq minimum number of sequences by OTUs by
#'   samples to take into count this OTUs in this sample. For example,
#'   if min_nb_seq=2,each value of 2 or less in the OTU table
#'   will not count in the venn diagram
#' @param taxa_fill (default NULL) fill the ASV upset using a column in
#'   `tax_table` slot.
#' @param na_remove : if TRUE (the default), NA values in fact are removed
#'   if FALSE, NA values are set to "NA"
#' @param numeric_fonction (default : sum) the function for numeric vector
#'   useful only for complex plot (see examples)
#' @param rarefy_after_merging Rarefy each sample after merging by the
#'   modalities of `fact` parameter
#' @param ... Other arguments passed on to the [ComplexUpset::upset()]
#'
#' @return A \code{\link{ggplot}}2 plot
#' @export
#' @author Adrien Taudière
#'
#' @seealso [ggvenn_pq()]
#' @examples
#' if (requireNamespace("ComplexUpset")) {
#'   upset_pq(data_fungi_mini,
#'     fact = "Height", width_ratio = 0.2,
#'     taxa_fill = "Class"
#'   )
#' }
#' \donttest{
#' if (requireNamespace("ComplexUpset")) {
#'   upset_pq(data_fungi_mini, fact = "Height", min_nb_seq = 1000)
#'   upset_pq(data_fungi_mini, fact = "Height", na_remove = FALSE)
#'
#'   upset_pq(data_fungi_mini, fact = "Time", width_ratio = 0.2, rarefy_after_merging = TRUE)
#'
#'   upset_pq(
#'     data_fungi_mini,
#'     fact = "Time",
#'     width_ratio = 0.2,
#'     annotations = list(
#'       "Sequences per ASV \n (log10)" = (
#'         ggplot(mapping = aes(y = log10(Abundance)))
#'         +
#'           geom_jitter(aes(
#'             color =
#'               Abundance
#'           ), na.rm = TRUE)
#'           +
#'           geom_violin(alpha = 0.5, na.rm = TRUE) +
#'           theme(legend.key.size = unit(0.2, "cm")) +
#'           theme(axis.text = element_text(size = 12))
#'       ),
#'       "ASV per phylum" = (
#'         ggplot(mapping = aes(fill = Phylum))
#'         +
#'           geom_bar() +
#'           ylab("ASV per phylum") +
#'           theme(legend.key.size = unit(0.2, "cm")) +
#'           theme(axis.text = element_text(size = 12))
#'       )
#'     )
#'   )
#'
#'   upset_pq(
#'     data_fungi_mini,
#'     fact = "Time",
#'     width_ratio = 0.2,
#'     numeric_fonction = mean,
#'     annotations = list(
#'       "Sequences per ASV \n (log10)" = (
#'         ggplot(mapping = aes(y = log10(Abundance)))
#'         +
#'           geom_jitter(aes(
#'             color =
#'               Abundance
#'           ), na.rm = TRUE)
#'           +
#'           geom_violin(alpha = 0.5, na.rm = TRUE) +
#'           theme(legend.key.size = unit(0.2, "cm")) +
#'           theme(axis.text = element_text(size = 12))
#'       ),
#'       "ASV per phylum" = (
#'         ggplot(mapping = aes(fill = Phylum))
#'         +
#'           geom_bar() +
#'           ylab("ASV per phylum") +
#'           theme(legend.key.size = unit(0.2, "cm")) +
#'           theme(axis.text = element_text(size = 12))
#'       )
#'     )
#'   )
#'
#'   upset_pq(
#'     subset_taxa(data_fungi_mini, Phylum == "Basidiomycota"),
#'     fact = "Time",
#'     width_ratio = 0.2,
#'     base_annotations = list(),
#'     annotations = list(
#'       "Sequences per ASV \n (log10)" = (
#'         ggplot(mapping = aes(y = log10(Abundance)))
#'         +
#'           geom_jitter(aes(
#'             color =
#'               Abundance
#'           ), na.rm = TRUE)
#'           +
#'           geom_violin(alpha = 0.5, na.rm = TRUE) +
#'           theme(legend.key.size = unit(0.2, "cm")) +
#'           theme(axis.text = element_text(size = 12))
#'       ),
#'       "ASV per phylum" = (
#'         ggplot(mapping = aes(fill = Class))
#'         +
#'           geom_bar() +
#'           ylab("ASV per Class") +
#'           theme(legend.key.size = unit(0.2, "cm")) +
#'           theme(axis.text = element_text(size = 12))
#'       )
#'     )
#'   )
#'
#'   data_fungi2 <- data_fungi_mini
#'   data_fungi2@sam_data[["Time_0"]] <- data_fungi2@sam_data$Time == 0
#'   data_fungi2@sam_data[["Height__Time_0"]] <-
#'     paste0(data_fungi2@sam_data[["Height"]], "__", data_fungi2@sam_data[["Time_0"]])
#'   data_fungi2@sam_data[["Height__Time_0"]][grepl("NA", data_fungi2@sam_data[["Height__Time_0"]])] <-
#'     NA
#'   upset_pq(data_fungi2, fact = "Height__Time_0", width_ratio = 0.2, min_size = 2)
#' }
#' }
upset_pq <- function(physeq,
                     fact,
                     taxa_fill = NULL,
                     min_nb_seq = 0,
                     na_remove = TRUE,
                     numeric_fonction = sum,
                     rarefy_after_merging = FALSE,
                     ...) {
  if (!is.null(min_nb_seq)) {
    physeq <- subset_taxa_pq(physeq, taxa_sums(physeq) >= min_nb_seq)
  }

  if (na_remove) {
    physeq <-
      subset_samples_pq(physeq, !is.na(physeq@sam_data[[fact]]))
  } else {
    physeq@sam_data[[fact]][is.na(physeq@sam_data[[fact]])] <-
      "NA"
  }

  physeq <- merge_samples2(physeq, fact)

  if (rarefy_after_merging) {
    physeq <- clean_pq(rarefy_even_depth(physeq))
  }

  psm <- psmelt(physeq)
  samp_names <- unique(psm$Sample)
  psm <-
    psm %>%
    mutate(val = TRUE) %>%
    tidyr::pivot_wider(names_from = Sample, values_from = val)
  psm[samp_names][is.na(psm[samp_names])] <- FALSE

  psm <- psm %>% filter(Abundance != 0)
  psm[[fact]] <- as.character(psm[[fact]])

  psm2 <- data.frame(lapply(psm, function(col) {
    tapply(col, paste0(psm$OTU), function(vec) {
      diff_fct_diff_class(vec, numeric_fonction = numeric_fonction, na.rm = TRUE)
    })
  })) %>% arrange(., desc(Abundance))

  colnames(psm2) <- colnames(psm)

  if (is.null(taxa_fill)) {
    p <-
      ComplexUpset::upset(psm2, intersect = samp_names, ...) + xlab(fact)
  } else {
    p <- ComplexUpset::upset(
      psm2,
      intersect = samp_names,
      base_annotations = list(),
      annotations = list("Taxa" = (
        ggplot(mapping = aes(fill = .data[[taxa_fill]])) +
          geom_bar() +
          ylab("Taxa per Class") +
          theme(legend.key.size = unit(0.2, "cm")) +
          theme(axis.text = element_text(size = 12))
      )),
      ...
    ) + xlab(fact)
  }

  return(p)
}
################################################################################


################################################################################
#' Test for differences between intersections
#'
#' @description 
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#' 
#' See [upset_pq()] to plot upset.
#'
#' @inheritParams upset_pq
#' @param var_to_test (default c("OTU")) : a vector of column present in
#'   the tax_table slot from the physeq object
#' @param ... Other arguments passed on to the [ComplexUpset::upset_test()]
#'
#' @return A \code{\link{ggplot}}2 plot
#' @export
#' @author Adrien Taudière
#'
#' @seealso [upset_pq()]
#' @examples
#' data(data_fungi)
#' if (requireNamespace("ComplexUpset")) {
#'   upset_test_pq(data_fungi, "Height", var_to_test = c("OTU", "Class", "Guild"))
#'   upset_test_pq(data_fungi, "Time")
#' }
upset_test_pq <-
  function(physeq,
           fact,
           var_to_test = "OTU",
           min_nb_seq = 0,
           na_remove = TRUE,
           numeric_fonction = sum,
           ...) {
    if (!is.null(min_nb_seq)) {
      physeq <- subset_taxa_pq(physeq, taxa_sums(physeq) >= min_nb_seq)
    }

    if (na_remove) {
      physeq <-
        subset_samples_pq(physeq, !is.na(physeq@sam_data[[fact]]))
    } else {
      physeq@sam_data[[fact]][is.na(physeq@sam_data[[fact]])] <-
        "NA"
    }

    physeq <- merge_samples2(physeq, fact)

    psm <- psmelt(physeq)
    samp_names <- unique(psm$Sample)
    psm <-
      psm %>%
      mutate(val = TRUE) %>%
      tidyr::pivot_wider(names_from = Sample, values_from = val)
    psm[samp_names][is.na(psm[samp_names])] <- FALSE

    psm <- psm %>% filter(Abundance != 0)
    psm[[fact]] <- as.character(psm[[fact]])

    psm2 <- data.frame(lapply(psm, function(col) {
      tapply(col, paste0(psm$OTU), function(vec) {
        diff_fct_diff_class(vec,
          numeric_fonction = numeric_fonction,
          na.rm = TRUE
        )
      })
    })) %>% arrange(., desc(Abundance))

    colnames(psm2) <- colnames(psm)

    res_test <-
      ComplexUpset::upset_test(psm2[, c(var_to_test, samp_names)],
        intersect = samp_names, ...
      )

    return(res_test)
  }
################################################################################

################################################################################
#' Compute different functions for different class of vector.
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Mainly an internal function useful in "sapply(..., tapply)" methods
#'
#' @param x : a vector
#' @param numeric_fonction : a function for numeric vector. For ex. `sum` or `mean`
#' @param logical_method : A method for logical vector. One of :
#'   - TRUE_if_one (default)
#'   - NA_if_not_all_TRUE
#'   - FALSE_if_not_all_TRUE
#' @param character_method : A method for character vector (and factor). One of :
#'   - unique_or_na (default)
#'   - more_frequent
#'   - more_frequent_without_equality
#' @param ... Other arguments passed on to the numeric function (ex. na.rm=TRUE)
#' @return a single value
#' @export
#'
#' @examples
#'
#' diff_fct_diff_class(
#'   data_fungi@sam_data$Sample_id,
#'   numeric_fonction = sum,
#'   na.rm = TRUE
#' )
#' diff_fct_diff_class(
#'   data_fungi@sam_data$Time,
#'   numeric_fonction = mean,
#'   na.rm = TRUE
#' )
#' diff_fct_diff_class(
#'   data_fungi@sam_data$Height == "Low",
#'   logical_method = "TRUE_if_one"
#' )
#' diff_fct_diff_class(
#'   data_fungi@sam_data$Height == "Low",
#'   logical_method = "NA_if_not_all_TRUE"
#' )
#' diff_fct_diff_class(
#'   data_fungi@sam_data$Height == "Low",
#'   logical_method = "FALSE_if_not_all_TRUE"
#' )
#' diff_fct_diff_class(
#'   data_fungi@sam_data$Height,
#'   character_method = "unique_or_na"
#' )
#' diff_fct_diff_class(
#'   c("IE", "IE"),
#'   character_method = "unique_or_na"
#' )
#' diff_fct_diff_class(
#'   c("IE", "IE", "TE", "TE"),
#'   character_method = "more_frequent"
#' )
#' diff_fct_diff_class(
#'   c("IE", "IE", "TE", "TE"),
#'   character_method = "more_frequent_without_equality"
#' )
#' @author Adrien Taudière
diff_fct_diff_class <-
  function(x,
           numeric_fonction = mean,
           logical_method = "TRUE_if_one",
           character_method = "unique_or_na",
           ...) {
    if (is.character(x) || is.factor(x)) {
      if (length(unique(x)) == 1) {
        return(unique(x))
      } else if (character_method == "unique_or_na") {
        return(NA)
      } else if (character_method == "more_frequent") {
        return(names(sort(table(x), decreasing = TRUE)[1]))
      } else if (character_method == "more_frequent_without_equality") {
        if (sort(table(x), decreasing = TRUE)[1] == sort(table(x), decreasing = TRUE)[2]) {
          return(NA)
        } else {
          return(names(sort(table(x), decreasing = TRUE)[1]))
        }
      } else {
        stop(paste0(
          character_method,
          " is not a valid method for character_method params."
        ))
      }
    } else if (is.numeric(x)) {
      return(numeric_fonction(x, ...))
    } else if (is.logical(x)) {
      if (logical_method == "TRUE_if_one") {
        if (sum(x, na.rm = TRUE) > 0) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      }
      if (logical_method == "NA_if_not_all_TRUE") {
        if (sum(x, na.rm = TRUE) > 0 && sum(!x, na.rm = TRUE) == 0) {
          return(TRUE)
        } else if (sum(!x, na.rm = TRUE) > 0 &&
          sum(x, na.rm = TRUE) > 0) {
          return(NA)
        } else if (sum(!x, na.rm = TRUE) > 0 &&
          sum(x, na.rm = TRUE) == 0) {
          return(FALSE)
        }
      }
      if (logical_method == "FALSE_if_not_all_TRUE") {
        if (sum(x, na.rm = TRUE) > 0 && sum(!x, na.rm = TRUE) == 0) {
          return(TRUE)
        } else {
          return(FALSE)
        }
      } else {
        stop(paste0(
          logical_method,
          " is not a valid method for character_method params."
        ))
      }
    } else {
      stop("At least one column is neither numeric nor character or logical")
    }
  }
################################################################################


################################################################################
#' Plot the distribution of sequences or ASV in one taxonomic levels
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Graphical representation of distribution of taxonomy, optionnaly across a factor.
#' 
#' @inheritParams clean_pq
#' @param fact Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}.
#' @param taxa (default: 'Order') Name of the taxonomic rank of interest
#' @param percent_bar (default FALSE) If TRUE, the stacked bar fill all
#'   the space between 0 and 1. It just set position = "fill" in the
#'   `ggplot2::geom_bar()` function
#' @param nb_seq (logical; default TRUE) If set to FALSE, only the number of ASV
#'   is count. Concretely, physeq otu_table is transformed in a binary
#'   otu_table (each value different from zero is set to one)
#' @return A \code{\link{ggplot}}2 plot  with bar representing the number of sequence en each
#'   taxonomic groups
#' @export
#'
#' @examples
#'
#' data_fungi_ab <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)
#' tax_bar_pq(data_fungi_ab) + theme(legend.position = "none")
#' \donttest{
#' tax_bar_pq(data_fungi_ab, taxa = "Class")
#' tax_bar_pq(data_fungi_ab, taxa = "Class", percent_bar = TRUE)
#' tax_bar_pq(data_fungi_ab, taxa = "Class", fact = "Time")
#' }
#' @author Adrien Taudière
#' @seealso [plot_tax_pq()] and [multitax_bar_pq()]
#'
tax_bar_pq <-
  function(physeq,
           fact = "Sample",
           taxa = "Order",
           percent_bar = FALSE,
           nb_seq = TRUE) {
    if (!nb_seq) {
      physeq <- as_binary_otu_table(physeq)
    }
    psm <- psmelt(physeq)
    if (percent_bar) {
      ggplot(psm) +
        geom_bar(aes(x = .data[[fact]], fill = .data[[taxa]], y = Abundance),
          stat = "identity",
          position = "fill"
        )
    } else {
      ggplot(psm) +
        geom_bar(aes(x = .data[[fact]], fill = .data[[taxa]], y = Abundance),
          stat = "identity"
        )
    }
  }
################################################################################

################################################################################
################################################################################
#' Ridge plot of a phyloseq object
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Graphical representation of distribution of taxa across a factor using ridges. 
#' 
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot
#'    different lines
#' @param nb_seq (logical; default TRUE) If set to FALSE, only the number of ASV
#'   is count. Concretely, physeq `otu_table` is transformed in a binary
#'   `otu_table` (each value different from zero is set to one)
#' @param log10trans (logical, default TRUE) If TRUE,
#'   the number of sequences (or ASV if nb_seq = FALSE) is log10
#'   transformed.
#' @param tax_level The taxonomic level to fill ridges
#' @param ... Other params passed on to [ggridges::geom_density_ridges()]
#'
#' @return A \code{\link{ggplot}}2 plot  with bar representing the number of sequence en each
#'   taxonomic groups
#' @export
#' @author Adrien Taudière
#' @examples
#' if (requireNamespace("ggridges")) {
#'   ridges_pq(data_fungi_mini, "Time", alpha = 0.5, log10trans = FALSE) + xlim(c(0, 1000))
#' }
#' \donttest{
#' if (requireNamespace("ggridges")) {
#'   ridges_pq(data_fungi_mini, "Time", alpha = 0.5, scale = 0.9)
#'   ridges_pq(data_fungi_mini, "Sample_names", log10trans = TRUE) + facet_wrap("~Height")
#'
#'   ridges_pq(data_fungi_mini,
#'     "Time",
#'     jittered_points = TRUE,
#'     position = ggridges::position_points_jitter(width = 0.05, height = 0),
#'     point_shape = "|", point_size = 3, point_alpha = 1, alpha = 0.7,
#'     scale = 0.8
#'   )
#' }
#' }
ridges_pq <- function(physeq,
                      fact,
                      nb_seq = TRUE,
                      log10trans = TRUE,
                      tax_level = "Class",
                      ...) {
  psm <- psmelt(physeq)
  psm <- psm %>% filter(Abundance > 0)

  if (log10trans) {
    psm$Abundance <- log10(psm$Abundance)
  }
  if (nb_seq) {
    p <- ggplot(psm, aes(y = factor(.data[[fact]]), x = Abundance)) +
      ggridges::geom_density_ridges(aes(fill = .data[[tax_level]]), ...) +
      xlim(c(0, NA))
  } else {
    psm_asv <-
      psm %>%
      group_by(.data[[fact]], OTU, .data[[tax_level]]) %>%
      summarise("count" = n())

    p <-
      ggplot(psm_asv, aes(y = factor(.data[[fact]]), x = count)) +
      ggridges::geom_density_ridges(
        aes(fill = .data[[tax_level]]),
        ...
      ) +
      xlim(c(0, NA))
  }
  return(p)
}
################################################################################

################################################################################

################################################################################
#' Plot treemap of 2 taxonomic levels
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Note that lvl2need to be nested in lvl1
#'
#' @inheritParams clean_pq
#' @param lvl1 (required) Name of the first (higher) taxonomic rank of interest
#' @param lvl2 (required) Name of the second (lower) taxonomic rank of interest
#' @param nb_seq (logical; default TRUE) If set to FALSE, only the number of ASV
#'   is count. Concretely, physeq otu_table is transformed in a binary
#'   otu_table (each value different from zero is set to one)
#' @param log10trans (logical, default TRUE) If TRUE,
#'   the number of sequences (or ASV if nb_seq = FALSE) is log10
#'   transformed.
#' @param plot_legend (logical, default FALSE) If TRUE, plot che
#'   legend of color for lvl 1
#' @param ... Other arguments passed on to [treemapify::geom_treemap()] function.
#'
#' @return A ggplot2 graphic
#' @export
#'
#' @author Adrien Taudière
#' @examples
#' data(data_fungi_sp_known)
#' if (requireNamespace("treemapify")) {
#'   treemap_pq(
#'     clean_pq(subset_taxa(
#'       data_fungi_sp_known,
#'       Phylum == "Basidiomycota"
#'     )),
#'     "Order", "Class",
#'     plot_legend = TRUE
#'   )
#' }
#' \donttest{
#' if (requireNamespace("treemapify")) {
#'   treemap_pq(
#'     clean_pq(subset_taxa(
#'       data_fungi_sp_known,
#'       Phylum == "Basidiomycota"
#'     )),
#'     "Order", "Class",
#'     log10trans = FALSE
#'   )
#'   treemap_pq(
#'     clean_pq(subset_taxa(
#'       data_fungi_sp_known,
#'       Phylum == "Basidiomycota"
#'     )),
#'     "Order", "Class",
#'     nb_seq = FALSE, log10trans = FALSE
#'   )
#' }
#' }
treemap_pq <- function(physeq,
                       lvl1,
                       lvl2,
                       nb_seq = TRUE,
                       log10trans = TRUE,
                       plot_legend = FALSE,
                       ...) {
  if (!nb_seq) {
    physeq <- as_binary_otu_table(physeq)
  }

  psm <- psmelt(physeq) %>%
    filter(!is.na(.data[[lvl2]])) %>%
    filter(!is.na(.data[[lvl1]]))

  psm2 <- psm %>%
    group_by(.data[[lvl2]]) %>%
    reframe(Abundance = sum(Abundance), LVL1 = unique(.data[[lvl1]]))

  if (log10trans) {
    psm2$Abundance <- log10(psm2$Abundance)
  }

  p <-
    ggplot(
      psm2,
      aes(
        area = Abundance,
        fill = LVL1,
        label = .data[[lvl2]],
        subgroup = LVL1
      )
    ) +
    treemapify::geom_treemap(...) +
    treemapify::geom_treemap_subgroup_border(colour = "white", size = 4) +
    treemapify::geom_treemap_text(
      colour = "white",
      place = "centre",
      size = 15,
      grow = TRUE
    )

  if (!plot_legend) {
    p <- p + theme(legend.position = "none")
  }

  if (nb_seq) {
    if (log10trans) {
      p <-
        p + ggtitle(paste0(
          "Nb of sequences (log10 transformed) by ",
          lvl1,
          " and ",
          lvl2
        ))
    } else {
      p <- p + ggtitle(paste0("Nb of sequences by ", lvl1, " and ", lvl2))
    }
  } else {
    if (log10trans) {
      p <- p + ggtitle(paste0(
        "Nb of ASV (log10 transformed) by ",
        lvl1, " and ", lvl2
      ))
    } else {
      p <- p + ggtitle(paste0("Nb of ASV by ", lvl1, " and ", lvl2))
    }
  }

  return(p)
}
################################################################################


################################################################################
#' Plot the partition the variation of a phyloseq object
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Graphical representation of the partition of variation obtain with [var_par_pq()].
#' @param res_varpart (required) the result of the functions [var_par_pq()]
#'   or [var_par_rarperm_pq()]
#' @param cutoff The values below cutoff will not be displayed.
#' @param digits The number of significant digits.
#' @param digits_quantile The number of significant digits for quantile.
#' @param fill_bg Fill colours of ellipses.
#' @param show_quantiles Do quantiles are printed ?
#' @param filter_quantile_zero Do we filter out value with quantile encompassing
#'   the zero value?
#' @param show_dbrda_signif Do dbrda significance for each component is printed
#'   using *?
#' @param show_dbrda_signif_pval (float, `[0:1]`) The value under which the
#'  dbrda is considered significant.
#' @param alpha (int, `[0:255]`) Transparency of the fill colour.
#' @param id.size A numerical value giving the character expansion factor for the names of circles or ellipses.
#' @param min_prop_pval_signif_dbrda (float, `[0:1]`) Only used if using the
#'   result of [var_par_rarperm_pq()] function. The * for dbrda_signif is only add if
#'   at least `min_prop_pval_signif_dbrda` of permutations show significance.
#'
#' @return A plot
#' @export
#' @author Adrien Taudière
#' @seealso [var_par_rarperm_pq()], [var_par_pq()]
#' @examples
#' \donttest{
#' if (requireNamespace("vegan")) {
#'   data_fungi_woNA <- subset_samples(data_fungi, !is.na(Time) & !is.na(Height))
#'   res_var_9 <- var_par_rarperm_pq(
#'     data_fungi_woNA,
#'     list_component = list(
#'       "Time" = c("Time"),
#'       "Size" = c("Height", "Diameter")
#'     ),
#'     nperm = 9,
#'     dbrda_computation = TRUE
#'   )
#'   res_var_2 <- var_par_rarperm_pq(
#'     data_fungi_woNA,
#'     list_component = list(
#'       "Time" = c("Time"),
#'       "Size" = c("Height", "Diameter")
#'     ),
#'     nperm = 2,
#'     dbrda_computation = TRUE
#'   )
#'   res_var0 <- var_par_pq(data_fungi_woNA,
#'     list_component = list(
#'       "Time" = c("Time"),
#'       "Size" = c("Height", "Diameter")
#'     ),
#'     dbrda_computation = TRUE
#'   )
#'   plot_var_part_pq(res_var0, digits_quantile = 2, show_dbrda_signif = TRUE)
#'   plot_var_part_pq(res_var_9,
#'     digits_quantile = 2, show_quantiles = TRUE,
#'     show_dbrda_signif = TRUE
#'   )
#'   plot_var_part_pq(
#'     res_var_2,
#'     digits = 5,
#'     digits_quantile = 2,
#'     cutoff = 0,
#'     show_quantiles = TRUE
#'   )
#' }
#' }
#' @importFrom stats anova as.formula quantile
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `vegan::varpart()` if you
#'   use this function.
plot_var_part_pq <-
  function(res_varpart,
           cutoff = 0,
           digits = 1,
           digits_quantile = 2,
           fill_bg = c("seagreen3", "mediumpurple", "blue", "orange"),
           show_quantiles = FALSE,
           filter_quantile_zero = TRUE,
           show_dbrda_signif = FALSE,
           show_dbrda_signif_pval = 0.05,
           alpha = 63,
           id.size = 1.2,
           min_prop_pval_signif_dbrda = 0.95) {
    if (show_dbrda_signif_pval > 1 | show_dbrda_signif_pval < 0) {
      stop("show_dbrda_signif_pval value must be within the range [0-1]")
    }
    if (min_prop_pval_signif_dbrda > 1 |
      min_prop_pval_signif_dbrda < 0) {
      stop("show_dbrda_signif_pval value must be within the range [0-1]")
    }
    x <- res_varpart$part
    vals <- x$indfract$Adj.R.square
    is.na(vals) <- vals < cutoff
    vals <- round(vals, digits + 1)
    labs_text <- format(vals, digits = digits, nsmall = digits + 1)
    labs_text <- gsub("NA", "", labs_text)
    if (show_quantiles) {
      labs_text <- paste0(
        labs_text,
        "\n (",
        round(x$indfract$Adj.R.squared_quantil_min, digits_quantile + 1),
        "...",
        round(x$indfract$Adj.R.squared_quantil_max, digits_quantile + 1),
        ")"
      )
      labs_text[is.na(vals)] <- ""
    }

    if (filter_quantile_zero) {
      labs_text[x$indfract$Adj.R.squared_quantil_min < 0] <- ""
    }

    if (show_dbrda_signif) {
      if (is.null(res_varpart$dbrda_result_prop_pval_signif)) {
        cond <-
          seq_along(res_varpart$dbrda_result)[sapply(res_varpart$dbrda_result, function(x) {
            x$`Pr(>F)`[[1]] < show_dbrda_signif_pval
          })]
        res_varpart$Xnames[cond] <-
          paste0(res_varpart$Xnames[cond], "*")
      } else {
        cond <-
          seq_along(res_varpart$dbrda_result)[res_varpart$dbrda_result_prop_pval_signif >=
            min_prop_pval_signif_dbrda]
        res_varpart$Xnames[cond] <-
          paste0(res_varpart$Xnames[cond], "*")
      }
    }

    vegan::showvarparts(
      x$nsets,
      labs_text,
      bg = fill_bg,
      alpha = alpha,
      id.size = id.size,
      Xnames = res_varpart$Xnames
    )
    if (any(is.na(vals))) {
      graphics::mtext(paste("Values <", cutoff, " not shown", sep = ""), 1)
    }
    if (sum(x$indfract$Adj.R.squared_quantil_min) > 0 &&
      filter_quantile_zero) {
      graphics::mtext(
        paste("Values with min quantile <0 not shown", sep = ""),
        side = 1,
        line = 1
      )
    }
    if (show_dbrda_signif) {
      if (is.null(res_varpart$dbrda_result_prop_pval_signif)) {
        graphics::mtext(
          paste(
            "* indicate significant anova of dbRDA for each component at p=",
            show_dbrda_signif_pval,
            sep = ""
          )
        )
      } else {
        graphics::mtext(
          paste(
            "* indicate significant anova of dbRDA, for each component, in at least ",
            round(min_prop_pval_signif_dbrda * 100, 2),
            "% of rarefaction permutations",
            sep = ""
          )
        )
      }
    }
    return(invisible())
  }
################################################################################


################################################################################
#' Scatterplot with marginal distributions and statistical results against
#' Hill diversity of phyloseq object
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Basically a wrapper of function [ggstatsplot::ggscatterstats()] for
#' object of class phyloseq and Hill number.
#'
#' @inheritParams clean_pq
#' @param num_modality (required) Name of the numeric column in
#'   `physeq@sam_data` to plot and test against hill numberk
#' @param hill_scales (a vector of integer) The list of q values to compute
#'   the hill number H^q. If Null, no hill number are computed. Default value
#'   compute the Hill number 0 (Species richness), the Hill number 1
#'   (exponential of Shannon Index) and the Hill number 2 (inverse of Simpson
#'   Index).
#' @param rarefy_by_sample (logical, default FALSE) If TRUE, rarefy
#'   samples using [phyloseq::rarefy_even_depth()] function.
#' @param one_plot (logical, default FALSE) If TRUE, return a unique
#'   plot with the three plot inside using the patchwork package.
#' @param ... Other arguments passed on to [ggstatsplot::ggscatterstats()]
#'   function.
#'
#' @return Either an unique ggplot2 object (if one_plot is TRUE) or
#'  a list of ggplot2 plot for each hill_scales.
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' if (requireNamespace("ggstatsplot")) {
#'   ggscatt_pq(data_fungi_mini, "Time", type = "non-parametric")
#'   ggscatt_pq(data_fungi_mini, "Time", hill_scales = 1:4, type = "parametric")
#'   ggscatt_pq(data_fungi_mini, "Sample_id",
#'     hill_scales = c(0, 0.5),
#'     one_plot = FALSE
#'   )
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `ggstatsplot::ggscatterstats()` if you
#'   use this function.
#' @seealso [ggbetween_pq()]
ggscatt_pq <- function(physeq,
                       num_modality,
                       hill_scales = c(0, 1, 2),
                       rarefy_by_sample = FALSE,
                       one_plot = TRUE,
                       ...) {
  verify_pq(physeq)
  physeq <- clean_pq(physeq, force_taxa_as_columns = TRUE)

  if (rarefy_by_sample) {
    physeq <- clean_pq(rarefy_even_depth(physeq))
  }

  p_list <- list()
  psm_res <- psmelt_samples_pq(physeq, hill_scales = hill_scales)
  for (i in seq_along(hill_scales)) {
    p_list[[i]] <-
      ggstatsplot::ggscatterstats(
        psm_res, !!paste0("Hill_", hill_scales[[i]]), !!num_modality,
        ...
      )
  }

  if (one_plot) {
    return(patchwork::wrap_plots(p_list))
  } else {
    return(p_list)
  }
}
################################################################################





################################################################################
#' Alluvial plot for taxonomy and samples factor vizualisation
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Basically a wrapper of [ggalluvial](https://corybrunson.github.io/ggalluvial/index.html)
#' package
#'
#' @inheritParams clean_pq
#' @param taxa_ranks A vector of taxonomic ranks. For examples c("Family","Genus").
#'   If taxa ranks is not set
#'   (default value = c("Phylum", "Class", "Order", "Family")).
#' @param wrap_factor A name to determine
#'   which samples to merge using [merge_samples2()] function.
#'   Need to be in \code{physeq@sam_data}.
#'   Need to be use when you want to wrap by factor the final plot
#'   with the number of taxa (type="nb_taxa")
#' @param by_sample (logical) If FALSE (default), sample information is not taking
#'   into account, so the taxonomy is studied globally. If fact is not NULL, by_sample
#'   is automatically set to TRUE.
#' @param rarefy_by_sample (logical, default FALSE) If TRUE, rarefy
#'   samples using [phyloseq::rarefy_even_depth()] function.
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot  the last column
#' @param type If "nb_seq" (default), the number of sequences is
#'   used in plot. If "nb_taxa", the number of ASV is plotted.
#' @param width (passed on to [ggalluvial::geom_flow()]) the width of each stratum,
#'   as a proportion of the distance between axes. Defaults to 1/3.
#' @param min.size (passed on to [ggfittext::geom_fit_text()]) Minimum font size,
#'   in points. Text that would need to be shrunk below this size to fit the box will
#'   be hidden. Defaults to 4 pt.
#' @param na_remove (logical, default FALSE) If set to TRUE, remove samples with
#'   NA in the variables set in formula.
#' @param use_ggfittext (logical, default FALSE) Do we use ggfittext to plot labels?
#' @param use_geom_label (logical, default FALSE) Do we use geom_label to plot labels?
#' @param size_lab Size for label if use_ggfittext is FALSE
#' @param ... Other arguments passed on to [ggalluvial::geom_flow()] function.
#'
#' @return A ggplot object
#' @export
#' @author Adrien Taudière
#' @examples
#' if (requireNamespace("ggalluvial")) {
#'   ggaluv_pq(data_fungi_mini)
#' }
#' \donttest{
#' if (requireNamespace("ggalluvial")) {
#'   ggaluv_pq(data_fungi_mini, type = "nb_taxa")
#'
#'   ggaluv_pq(data_fungi_mini, wrap_factor = "Height", by_sample = TRUE, type = "nb_taxa") +
#'     facet_wrap("Height")
#'
#'   ggaluv_pq(data_fungi_mini,
#'     width = 0.9, min.size = 10,
#'     type = "nb_taxa", taxa_ranks = c("Phylum", "Class", "Order", "Family", "Genus")
#'   ) +
#'     coord_flip() + scale_x_discrete(limits = rev)
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `ggalluvial` package if you
#'   use this function.
#' @seealso [sankey_pq()]
ggaluv_pq <- function(physeq,
                      taxa_ranks = c("Phylum", "Class", "Order", "Family"),
                      wrap_factor = NULL,
                      by_sample = FALSE,
                      rarefy_by_sample = FALSE,
                      fact = NULL,
                      type = "nb_seq",
                      width = 1.2,
                      min.size = 3,
                      na_remove = FALSE,
                      use_ggfittext = FALSE,
                      use_geom_label = FALSE,
                      size_lab = 2,
                      ...) {
  verify_pq(physeq)
  if (rarefy_by_sample) {
    physeq <- rarefy_even_depth(physeq)
  }

  if (na_remove && !is.null(fact)) {
    physeq <- subset_samples_pq(physeq, !is.na(physeq@sam_data[[fact]]))
  }

  if (!is.null(wrap_factor)) {
    physeq <-
      merge_samples2(physeq, physeq@sam_data[[wrap_factor]])
  } else if (!by_sample || !is.null(fact)) {
    physeq <-
      merge_samples2(physeq, group = rep("all_samples_together", nsamples(physeq)))
  }

  if (type == "nb_taxa") {
    physeq <- as_binary_otu_table(physeq)
  } else if (type != "nb_seq") {
    stop("Type must be eiter nb_seq or nb_taxa")
  }

  psm_samp <-
    psmelt_samples_pq(
      physeq,
      taxa_ranks = taxa_ranks,
      hill_scales = NULL,
      rarefy_by_sample = FALSE
    )

  if (is.null(fact)) {
    psm_samp <- ggalluvial::to_lodes_form(psm_samp, axes = taxa_ranks)
  } else {
    psm_samp <- ggalluvial::to_lodes_form(psm_samp, axes = c(taxa_ranks, fact))
  }

  p <- ggplot(
    data = psm_samp,
    aes(
      alluvium = alluvium,
      x = x,
      stratum = stratum,
      y = Abundance,
      fill = after_stat(stratum)
    )
  ) +
    ggalluvial::geom_flow(...) +
    ggalluvial::geom_stratum() +
    theme_minimal() +
    theme(legend.position = "none")

  if (use_ggfittext) {
    p <- p +
      ggfittext::geom_fit_text(
        aes(label = stratum),
        stat = "stratum",
        width = width,
        min.size = min.size
      )
  } else if (use_geom_label) {
    p <- p +
      geom_label(
        aes(label = stratum),
        stat = "stratum",
        size = size_lab
      )
  }

  if (!is.null(wrap_factor)) {
    p <- p + facet_wrap(wrap_factor)
  }
  return(p)
}
################################################################################
