################################################################################
#Plot the result of a mt test
################################################################################
#' Plot the result of a mt test (\code{\link[phyloseq]{mt}})
#' @description 
#' `r lifecycle::badge("maturing")`
#'  
#' @param mt (required): result of a mt test
#' @param alpha (default = 0.05): Choose the cut off p-value to plot taxa
#' @param color_tax : A taxonomic level to color the points
#' @param taxa : The taxonomic level choose for x-positioning
#' data("GlobalPatterns")
#' res <- mt(GlobalPatterns, "SampleType", test="f")
#' plot_mt(res, color_tax = "Phylum") + scale_color_hue()
#' @author Adrien Taudière
#'
#' @return a \code{\link{ggplot}}2 plot of result of a mt test
#' @export
#' @seealso \code{\link[phyloseq]{mt}}

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
      factor(d$tax, levels =
               unique(factor(as.character(d[, "Species"]))
                      [rev(order(d$teststat))]))

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
#' @description 
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (required): Name of the factor in physeq@sam_data used to plot
#'    different lines
#' @param nb_seq (logical): Either plot accumulation curves using sequences or
#'   using samples
#' @param step (integer): distance among points calculated to plot lines. A
#'  low value give better plot but is more time consuming.
#'  Only used if nb_seq = TRUE.
#' @param by.fact (logical): First merge the OTU table by factor to plot only
#'   one line by factor
#' @param ci_col : Color vector for confidence intervall.
#'   Only use if nb_seq = FALSE.
#'   If nb_seq = TRUE, you can use ggplot to modify the plot.
#' @param col : Color vector for lines. Only use if nb_seq = FALSE.
#'   If nb_seq = TRUE, you can use ggplot to modify the plot.
#' @param lwd  (default = 3): thickness for lines. Only use if nb_seq = FALSE.
#' @param leg (logical): Plot legend or not. Only use if nb_seq = FALSE.
#' @param print_sam_names (logical): Print samples names or not?
#'    Only use if nb_seq = TRUE.
#' @param ci (default = 2) : Confidence intervall value used to multiply the
#'   standard error to plot confidence intervall
#' @param ... Additional arguments passed on to \code{\link{ggplot}}
#' if nb_seq = TRUE or to \code{\link{plot}} if nb_seq = FALSE
#'
#' data("GlobalPatterns")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == 'Archaea')
#' accu_plot(GP, 'SampleType', nb_seq = TRUE,
#'  by.fact = TRUE)
#'
#' @return A \code{\link{ggplot}}2 plot representing the richness
#' accumulation plot if nb_seq = TRUE, else, if nb_seq = FALSE
#' return a base plot.
#'
#' @export
#' @author Adrien Taudière
#' @seealso \code{\link[vegan]{specaccum}}
accu_plot <-
  function(physeq,
           fact = NULL,
           nb_seq = TRUE,
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

    if (!nb_seq) {
      factor_interm <-
        eval(parse(text = paste("physeq@sam_data$", fact, sep = "")))
      factor_interm <- as.factor(factor_interm)

      physeq_accu <- as.matrix(t(physeq@otu_table))
      physeq_accu[physeq_accu > 0] <- 1
      accu_all <- vegan::specaccum(physeq_accu)

      accu <- list()
      for (i in 1:nlevels(factor_interm)) {
        accu[[i]] <-
          vegan::specaccum(physeq_accu[factor_interm ==
                                         levels(factor_interm)[i], ])
      }

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

      if (is.null(col)) {
        col <- funky_color(nlevels(factor_interm) + 1)
      }
      if (is.null(ci_col)) {
        transp <- function(col, alpha = 0.5) {
          res <-
            apply(grDevices::col2rgb(col), 2, function(c)
              grDevices::rgb(c[1] / 255, c[2] / 255, c[3] / 255, alpha))
          return(res)
        }
        ci_col <-
          transp(funky_color(nlevels(factor_interm) + 1), 0.3)
      }

      plot(
        accu_all,
        ci_type = "poly",
        ci_col = ci_col[1],
        col = col[1],
        lwd = lwd,
        ci_lty = 0,
        xlab = "Sample",
        ...
      )

      for (i in 1:nlevels(factor_interm)) {
        graphics::lines(
          accu[[i]],
          ci_type = "poly",
          ci_col = ci_col[i + 1],
          col = col[i + 1],
          lwd = lwd,
          ci_lty = 0
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

    if (nb_seq) {
      fact_interm <-
        as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))

      if (!by.fact) {
        x <- t(physeq@otu_table)
      } else {
        x <- apply(physeq@otu_table, 1, function(x)
          tapply(x, fact_interm, sum))
      }

      tot <- rowSums(x)
      nr <- nrow(x)

      if (is.null(step)) {
        step <- round(max(tot) / 30, 0)
      }
      if (is.null(step)) {
        step <- 1
      }

      n_max <- seq(1, max(tot), by = step)
      out <- lapply(seq_len(nr), function(i) {
        res <- list()
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) {
          n <- c(n, tot[i])
        }
        res_interm <- rarefy(x[i, ], n, se = TRUE)
        res <-
          cbind(as.matrix(res_interm)[1, ], as.matrix(res_interm)[2, ])
        return(res)
      })

      names(out) <- names(tot)

      df <- plyr::ldply(out, data.frame)

      cond <- c()
      for (i in 1:nlevels(as.factor(df$.id))) {
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
      dff$.id <-  names(dff$ylab)
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
        geom_line() + xlab("Number of sequences") +
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
#' Plot OTU circle for \code{\link{phyloseq-class}} object
#' @description 
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (required): Name of the factor to cluster samples by modalities.
#'        Need to be in \code{physeq@sam_data}.
#' @param taxa (Default:'Order'): Name of the taxonomic rank of interest
#' @param nb_seq (Default: TRUE): Represent the number of sequences or the
#'    number of OTUs (nb_seq = FALSE)
#' @param rarefy (logical): Does each samples modalities need to be rarefy in
#'               order to compare them with the same amount of sequences?
#' @param min_prop_tax (Default: 0.01): The minimum proportion for taxon to be
#'                     ploted
#' @param min_prop_mod (Default: 0.1) : The minimum proportion for modalities
#'                     to be ploted
#' @param gap_degree : Gap between two neighbour sectors.
#'                    It can be a single value or a vector. If it is a vector,
#'                     the first value corresponds to the gap after the first
#'                     sector.
#' @param start_degree : The starting degree from which the circle begins to
#'   draw. Note this degree is measured in the standard polar coordinate
#'   which means it is always reverse-clockwise.
#' @param row_col : Color vector for row
#' @param grid_col : Grid colors which correspond to sectors. The length of the
#'    vector should be either 1 or the number of sectors.
#'    It's preferred that grid_col is a named vector of which names
#'    correspond to sectors. If it is not a named vector, the
#'    order of grid_col corresponds to order of sectors.
#' @param log10trans (logical): Should sequence be log10 transformed
#'                   (more precisely by log10(1+x))?
#' @param ... Additional arguments passed on to
#'   \code{\link[circlize]{chordDiagram}} or \code{\link[circlize]{circos.par}}
#'
#' @examples
#' data("GlobalPatterns")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == 'Archaea')
#' otu_circle(GP, 'SampleType')
#' otu_circle(GP, 'SampleType', nb_seq = FALSE)
#' otu_circle(GP, 'SampleType', taxa = "Class")
#' @author Adrien Taudière
#'
#' @return A \code{\link[circlize]{chordDiagram}} plot representing the
#'   distribution of OTUs or sequences in the different modalities of the factor
#'   fact
#'
#' @export
#' @seealso \code{\link[circlize]{chordDiagram}}
#' @seealso \code{\link[circlize]{circos.par}}

otu_circle <-
  function(physeq = NULL,
           fact = NULL,
           taxa = "Order",
           nb_seq = TRUE,
           rarefy = FALSE,
           min_prop_tax = 0.01,
           min_prop_mod = 0.1,
           gap_degree = NULL,
           start_degree = NULL,
           row_col = NULL,
           grid_col = NULL,
           log10trans = F,
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }

    if (!physeq@otu_table@taxa_are_rows) {
      otu_tab <- t(physeq@otu_table)
    } else {
      otu_tab <- physeq@otu_table
    }

    if (!nb_seq) {
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
      apply(otu_tab, 2, function(x)
        tapply(x, physeq@tax_table[, taxcol],
               function(xx)
                 sum(xx, na.rm = T)))
    otu_table_ech <-
      apply(otu_table_tax, 1, function(x)
        tapply(x, physeq@sam_data[, taxsamp],
               function(xx)
                 sum(xx, na.rm = T)))
    if (rarefy) {
      otu_table_ech_interm <-
        rrarefy(otu_table_ech, min(rowSums(otu_table_ech)))
      message(
        paste(
          "Rarefaction by modalities deletes ",
          sum(otu_table_ech) - sum(otu_table_ech_interm),
          " (",
          round(100 * (
            sum(otu_table_ech) - sum(otu_table_ech_interm)
          ) / sum(otu_table_ech),
          2),
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
      otu_table_ech[(rowSums(otu_table_ech) / sum(otu_table_ech)) >
                      min_prop_mod,
                    (colSums(otu_table_ech) / sum(otu_table_ech)) >
                      min_prop_tax]
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
          round(100 * ncol(o_t_e_interm) / ncol(otu_table_ech),
                2),
          "%). Use 'min_prop_tax' to plot more taxa",
          sep = ""
        )
      )
    }
    otu_table_ech <- o_t_e_interm

    if (log10trans) {
      otu_table_ech <- apply(otu_table_ech, 2, function(x)
        log10(1 + x))
    }


    if (is.null(gap_degree)) {
      col2keep <- rep(1, ncol(otu_table_ech) - 1)
      row2keep <- rep(1, nrow(otu_table_ech) - 1)
      gap_degree <- c(row2keep, 10, col2keep, 10)

    }
    if (is.null(start_degree)) {
      start_degree <- 170
    }

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

    if (is.null(grid_col)) {
      grid_col <-
        c(funky_color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
    }

    if (is.null(row_col)) {
      row_col <-
        c(funky_color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
    }

    circlize::circos.par(gap.degree = gap_degree, start.degree = start_degree,
                         ...)
    circlize::chordDiagram(otu_table_ech, row.col = row_col,
                           grid.col = grid_col, ...)
    circlize::circos.clear()
  }
################################################################################

################################################################################
#' Sankey plot of \code{\link{phyloseq-class}} object
#' @description 
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (Optional): Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}.
#' @param taxa (Default: c(1:4)): a vector of taxonomic rank to plot
#' @param nb_seq (Default: FALSE): Represent the number of sequences or the
#'   number of OTUs (nb_seq = FALSE). Note that ploting the number of sequences
#'   is slower.
#' @param min_prop_tax (Default: 0): The minimum proportion for taxon to be
#'  ploted. EXPERIMENTAL. For the moment each links below the min.prop.
#'  tax is discard from the sankey network resulting in sometimes weird plot.
#' @param tax2remove : a vector of taxonomic groups to remove from the analysis
#'   (e.g. \code{c('Incertae sedis', 'unidentified')})
#' @param units : character string describing physical units (if any) for Value
#' @param symbol2sub (default = c('\\.', '-')): vector of symbol to delete in
#'   the taxonomy
#' @param ... Additional arguments passed on to
#'   \code{\link[networkD3]{sankeyNetwork}}
#'
#' @examples
#' data("GlobalPatterns")
#' GP <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[,1] == "Archaea")
#' sankey_phyloseq(GP, fact = 'SampleType')
#' sankey_phyloseq(GP, taxa = c(1:4), min_prop_tax = 0.01)
#' sankey_phyloseq(GP, taxa = c(1:4), min_prop_tax = 0.01, nb_seq = TRUE)
#' @author Adrien Taudière
#'
#' @return A \code{\link[networkD3]{sankeyNetwork}} plot representing the
#'  taxonomic distribution of OTUs or sequences. If \code{fact} is set,
#'  represent the distribution of the last taxonomic level in the modalities
#'  of \code{fact}
#'
#' @export
#' @seealso \code{\link[networkD3]{sankeyNetwork}}

sankey_phyloseq <-
  function(physeq = NULL,
           fact = NULL,
           taxa = c(1:4),
           nb_seq = FALSE,
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

    if (!nb_seq) {
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
    } else if (nb_seq) {
      mat_interm <- matrix()
      mat <- matrix(ncol = 3)
      colnames(mat) <- c("Var1", "Var2", "value")
      tax_table_interm <-
        physeq@tax_table[rep(1:ntaxa(physeq),
                             times = taxa_sums(physeq))]

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
        for (i in 1:dim(m)[1]) {
          for (j in 1:dim(m)[2]) {
            if (m[i, j] > 0) {
              res <- rbind(res, c(rownames(m)[i], colnames(m)[j], m[i, j]))
            }
          }
        }
        return(res)
      }

      mat_interm <-
        apply(otu_tab, 1, function(x)
          tapply(x, physeq@sam_data[, fact],
                 sum))

      if (!nb_seq) {
        mat_interm <-
          apply(mat_interm, 1, function(x)
            tapply(x, physeq@tax_table[,
                                       taxa[length(taxa)]], function(x)
                                         sum(x > 0)))
      } else if (nb_seq) {
        mat_interm <-
          apply(mat_interm, 1, function(x)
            tapply(x, physeq@tax_table[,
                                       taxa[length(taxa)]], sum))
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
      mat <- apply(mat, 2, function(x)
        gsub(symbol2sub[i], "", x))
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
        gsub(paste("\\<", tax_sank$nodes[i, 2], "\\>", sep = ""),
             tax_sank$nodes[i,
                            1],
             mat2[, 1])
      mat2[, 2] <-
        gsub(paste("\\<", tax_sank$nodes[i, 2], "\\>", sep = ""),
             tax_sank$nodes[i,
                            1],
             mat2[, 2])
    }

    tax_sank$links <- apply(mat2, 2, as.numeric)
    tax_sank$links <-
      data.frame(tax_sank$links[rowSums(is.na(tax_sank$links)) == 0, ])
    tax_sank$nodes <-
      as.data.frame(as.character(tax_sank$nodes[, 2]))
    names(tax_sank$nodes) <- c("name")
    names(tax_sank$links) <- c("source", "target", "value")
    if (is.null(units)) {
      if (!nb_seq) {
        units <- "OTUs"
      } else if (nb_seq) {
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
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (required): Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}.
#' @param min_nb_seq (Default: 0)): minimum number of sequences by OTUs by
#'  samples to take into count this OTUs in this sample. For example,
#'  if min_nb_seq=2,each value of 2 or less in the OTU table
#'  will be change into 0 for the analysis
#' @param print_values (logical) : Print (or not) the table of number of OTUs
#' for each combination.
#' If print_values is TRUE the object is not a ggplot object.
#' Please use print_values = FALSE if you want to add ggplot function
#' (cf example).
#'
#'
#' data("enterotype")
#' venn_phyloseq(enterotype, fact = 'SeqTech')
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus')
#' venn_phyloseq(enterotype, fact = 'Nationality', print_values = F)
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus', print_values = F) +
#' scale_fill_hue()
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus', print_values = F) +
#' scale_fill_hue()
#'
#' @return A \code{\link{ggplot}}2 plot representing Venn diagramm of
#' modalities of the argument \code{factor}
#'
#' @export
#' @author Adrien Taudière
#' @seealso \code{\link[venneuler]{venneuler}}

venn_phyloseq <-
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
        t(apply(physeq@otu_table, 1, function(x)
          by(x, moda, max)))
    } else if (length(moda) != dim(physeq@otu_table)[2]) {
      data_venn <-
        t(apply(t(physeq@otu_table), 1, function(x)
          by(x, moda, max)))
    } else {
      stop("The factor length and the number of samples must be identical")
    }
    combinations <- data_venn > min_nb_seq

    e <- new.env(TRUE, emptyenv())
    cn <- colnames(combinations)
    for (i in seq.int(dim(combinations)[1]))
      if (any(combinations[i, ])) {
        ec <- paste(cn[combinations[i, ]], collapse = "&")
        e[[ec]] <- if (is.null(e[[ec]]))
          1L
        else
          (e[[ec]] + 1L)
      }

    en <- ls(e, all.names = TRUE)
    weights <- as.numeric(unlist(lapply(en, get, e)))
    combinations <- as.character(en)

    table_value <-
      data.frame(combinations = as.character(combinations),
                 weights = as.double(weights))

    venn <- venneuler::venneuler(data_venn > min_nb_seq)
    venn_res <-
      data.frame(
        x = venn$centers[, 1],
        y = venn$centers[, 2],
        radius = venn$diameters / 2
      )

    nmod <- nrow(venn_res)
    x1 <- list()
    for (i in 1:nmod) {
      x1[[i]] <- grep(rownames(venn_res)[i], table_value$combinations)
    }

    for (i in seq_len(nrow(table_value))) {
      table_value$x[i] <-
        mean(venn$centers[, "x"][unlist(lapply(x1,
                                               function(x)
                                                 sum(x %in% i) > 0))])
      table_value$y[i] <-
        mean(venn$centers[, "y"][unlist(lapply(x1,
                                               function(x)
                                                 sum(x %in% i) > 0))])
    }

    df <- venn_res
    df$xlab <- df$x + (df$x - mean(df$x))
    df$ylab <- df$y + (df$y - mean(df$y))

    circularise <- function(d, n = 360) {
      angle <- seq(-pi, pi, length = n)
      make_circle <- function(x, y, r, modality) {
        data.frame(x = x + r * cos(angle),
                   y = y + r * sin(angle),
                   modality)
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
      ggplot() + geom_polygon(data = circles,
                              aes(x, y, group = modality, fill = modality),
                              alpha = 0.5) + theme_void()

    if (print_values) {
      g_legend <- function(agplot) {
        tmp <- ggplot_gtable(ggplot_build(agplot))
        leg <-
          which(sapply(tmp$grobs, function(x)
            x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
      }
      legend <- g_legend(p)


      grid::grid.newpage()
      vp1 <- viewport(
        width = 0.75,
        height = 1,
        x = 0.375,
        y = .5
      )
      vpleg <-
        viewport(
          width = 0.25,
          height = 0.5,
          x = 0.85,
          y = 0.75
        )
      subvp <- viewport(
        width = 0.3,
        height = 0.3,
        x = 0.85,
        y = 0.25
      )
      print(p + theme(legend.position = "none"), vp = vp1)
      upViewport(0)
      pushViewport(vpleg)
      grid.draw(legend)
      #Make the new viewport active and draw
      upViewport(0)
      pushViewport(subvp)
      grid.draw(gridExtra::tableGrob(table_value[, c(1, 2)], rows = NULL))
    }
    else{
      return(p)
    }
  }
################################################################################

################################################################################
#' Multiple plot function
#' @description 
#' `r lifecycle::badge("stable")`
#'  
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot
# objects)
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
#' @param ... : list of ggplot objects
#' @param plotlist : list of ggplot objects
#' @param cols : number of columns
#' @param layout : A matrix specifying the layout. If present, 'cols' is ignored.

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
                       nrow = ceiling(num_plots / cols))
    }

    if (num_plots == 1) {
      message(plots[[1]])

    } else {
      # Set up the page
      grid::grid.newpage()
      grid::pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                       ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:num_plots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <-
          as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]],
              vp = viewport(
                layout.pos.row = matchidx$row,
                layout.pos.col = matchidx$col
              ))
      }
    }
  }
################################################################################

################################################################################
#' Graphical representation of hill number 0, 1 and 2 accross a factor
#' @description 
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): A \code{\link{phyloseq-class}} object
#' @param variable (required): The variable to test
#' @param color_fac (optional): The variable to color the barplot
#' @param letters (optional, default=FALSE): If set to TRUE, the plot
#' show letters based on p-values for comparison. Use the
#'  \code{\link[multcompView]{multcompLetters}} function from the package
#'  multcompLetters.
#'
#' @return A list of 4 ggplot2 plot.
#' - plot_Hill_0 : the boxplot of Hill number 0 (= species richness)
#'     against the variable
#' - plot_Hill_1 : the boxplot of Hill number 1 (= Shannon index)
#'      against the variable
#' - plot_Hill_2 : the boxplot of Hill number 2 (= Simpson index)
#'     against the variable
#' - plot_tuckey : plot the result of the Tuckey HSD test
#'
#' @export

hill_phyloseq <-
  function(physeq,
           variable,
           color_fac = NA,
           letters = FALSE) {
    var <- sym(variable)
    if (is.na(color_fac)) {
      color_fac <- sym(variable)
    }
    else{
      color_fac <- sym(color_fac)
    }

    otu_hill <-
      vegan::renyi(physeq@otu_table, scale = c(0, 1, 2), hill = T)
    colnames(otu_hill) <- c("Hill_0", "Hill_1", "Hill_2")

    df_hill <- data.frame(otu_hill, physeq@sam_data)
    df_hill[, c(1:3)] <- apply(df_hill[, c(1:3)], 2, as.numeric)


    p_var <- hill_tuckey_phyloseq(physeq, variable)

    p_0 <- ggplot(df_hill, aes(group = !!var, Hill_0))  +
      geom_boxplot(outlier.size = 2, aes(colour = as.factor(!!color_fac)))
    p_1 <- ggplot(df_hill, aes(group = !!var, Hill_1))  +
      geom_boxplot(outlier.size = 2, aes(colour = as.factor(!!color_fac)))
    p_2 <- ggplot(df_hill, aes(group = !!var, Hill_2))  +
      geom_boxplot(outlier.size = 2, aes(colour = as.factor(!!color_fac)))


    if (letters) {
      ### HILL 0
      data_h0 <-
        p_var$data[grep("Hill Number 0", p_var$data[, 5]),]
      data_h0_pval <- data_h0$p.adj
      names(data_h0_pval) <- data_h0$modality
      letters <- multcompLetters(data_h0_pval, reversed = T)$Letters

      p_0 <- p_0 +
        geom_label(
          data = p_0$data %>% group_by(!!var) %>%
            summarise(max_Hill = max(Hill_0)),
          aes(x = max_Hill + 1),
          label = letters[match(names(letters), levels(factor(as.matrix(
            physeq@sam_data[, variable]
          ))))],
          y = ggplot_build(p_0)$data[[1]]$y,
          size = 4,
          stat = "unique",
          parse = TRUE
        )

      ### HILL 1

      data_h1 <-
        p_var$data[grep("Hill Number 1", p_var$data[, 5]), ]
      data_h1_pval <- data_h1$p.adj
      names(data_h1_pval) <- data_h1$modality
      letters <- multcompLetters(data_h1_pval, reversed = T)$Letters

      p_1 <- p_1 +
        geom_label(
          data = p_1$data %>% group_by(!!var) %>%
            summarise(max_Hill = max(Hill_1)),
          aes(x = max_Hill + 1),
          label = letters[match(names(letters), levels(factor(as.matrix(
            physeq@sam_data[, variable]
          ))))],
          y = ggplot_build(p_1)$data[[1]]$y,
          size = 4,
          stat = "unique",
          parse = TRUE
        )

      ### HILL 2

      data_h2 <-
        p_var$data[grep("Hill Number 2", p_var$data[, 5]),]
      data_h2_pval <- data_h2$p.adj
      names(data_h2_pval) <- data_h2$modality
      letters <- multcompLetters(data_h2_pval, reversed = T)$Letters

      p_2 <- p_2 +
        geom_label(
          data = p_2$data %>% group_by(!!var) %>%
            summarise(max_Hill = max(Hill_2)),
          aes(x = max_Hill + 0.5),
          label = letters[match(names(letters), levels(factor(as.matrix(
            physeq@sam_data[, variable]
          ))))],
          y = ggplot_build(p_1)$data[[1]]$y,
          size = 4,
          stat = "unique",
          parse = TRUE
        )

    }

    res <- list(
      "plot_Hill_0" = p_0,
      "plot_Hill_1" = p_1,
      "plot_Hill_2" = p_2,
      "plot_tuckey" = p_var
    )
  }
################################################################################

################################################################################
#' Make a datatable with the taxonomy of a \code{\link{phyloseq-class}} object
#' @description 
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): A \code{\link{phyloseq-class}} object
#' @param abundance (Default: TRUE): Does the number of sequences is print
#' @param taxonomic_level (Default: NULL): a vector of selected taxonomic
#' level using their column numbers (e.g. taxonomic_level = c(1:7))
#' @param ... Other argument for the datatable function
#'
#' @author Adrien Taudière
#' @return A datatable
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' tax_datatable(subset_taxa(GlobalPatterns,
#'               rowSums(GlobalPatterns@otu_table)>10000))
tax_datatable <- function(physeq, abundance = TRUE, taxonomic_level=NULL, ...) {
  df <- as.data.frame(physeq@tax_table)

  if (!is.null(taxonomic_level)){
    df <- df[,taxonomic_level]
  }

  if (abundance) {
    if (physeq@otu_table@taxa_are_rows) {
      df$nb_seq <- rowSums(physeq@otu_table)
    } else {
      df$nb_seq <- colSums(physeq@otu_table)
    }
  }

  DT::datatable(df, ...) %>%  DT::formatStyle(
    "nb_seq",
    background = DT::styleColorBar(df$nb_seq, "steelblue"),
    backgroundSize = "100% 90%",
    backgroundRepeat = "no-repeat",
    backgroundPosition = "center"
  )
}
################################################################################

################################################################################
#' Summarise a \code{\link{phyloseq-class}} object using a plot.
#' @description 
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): A \code{\link{phyloseq-class}} object
#' @examples
#'   data(esophagus)
#'   summary_plot_phyloseq(esophagus)
#' @return A ggplot2 object
#' @export
summary_plot_phyloseq <- function(physeq) {
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
        "occurences"
      ),
      paste(ncol(physeq@tax_table), "taxonomic \n levels"),
      paste(
        ncol(physeq@sam_data),
        "variables\n",
        nsamples(physeq),
        "samples"
      ),
      paste("Sequences length:\n",
            round(mean(
              Biostrings::width(physeq@refseq)
            ), 2), "+/-",
            round(stats::sd(
              Biostrings::width(physeq@refseq)
            ), 2))
    )
  )
  ggplot() +
    scale_x_continuous(name = "x") +
    scale_y_reverse(name = "y")  +
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
    geom_text(data = d,
              aes(
                x = x1 + (x2 - x1) / 2,
                y = y1 + (y2 - y1) / 1.7,
                label = nb_values
              ),
              size = 5) +
    geom_text(data = d,
              aes(
                x = x1 + (x2 - x1) / 2,
                y = y1 + (y2 - y1) / 5,
                label = names
              ),
              size = 7) +
    annotate(
      geom = "text",
      x = 0.65,
      y = 3,
      label = "Taxa",
      size = 7,
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
      size = 7,
      fontface = 2,
      color = "khaki4"
    )
}
################################################################################


#' Heat tree from `metacoder` package using `tax_table` slot
#' @description 
#' `r lifecycle::badge("maturing")`
#' @param physeq (required): A \code{\link{phyloseq-class}} object
#' @param taxonomic_level (Default: NULL): a vector of selected taxonomic level using
#'   their column numbers (e.g. taxonomic_level = c(1:7))
#' @param ... : Arguments parsed to \code{\link[metacoder]{heat_tree}}
#'
#' @return A plot
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' data("GlobalPatterns")
#' GPsubset <- subset_taxa(GlobalPatterns,
#'                         GlobalPatterns@tax_table[, 1] == 'Bacteria')
#' GPsubset <- subset_taxa(GPsubset,
#'                         rowSums(GPsubset@otu_table) > 1000)
#' GPsubset <- subset_taxa(GPsubset, rowSums(is.na(GPsubset@tax_table)) == 0)
#'
#' physeq_heat_tree(GPsubset,
#'                  node_size = n_obs,
#'                  node_color = n_obs,
#'                  node_label = taxon_names,
#'                  tree_label = taxon_names,
#'                  node_size_trans = "log10 area")

physeq_heat_tree <- function(physeq, taxonomic_level=NULL, ...) {
  if (!is.null(taxonomic_level)) {
    physeq@tax_table <- physeq@tax_table[,taxonomic_level]
  }
  
  data_metacoder <- metacoder::parse_phyloseq(physeq)
  metacoder::heat_tree(data_metacoder, ...)
}
