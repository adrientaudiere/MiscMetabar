####################################################################################
#Plot the result of a mt test
####################################################################################
#' Plot the result of a mt test (\code{\link[phyloseq]{mt}})
#' @param mt (required): result of a mt test
#' @param alpha (default = 0.05): Choose the cut off p-value to plot taxa
#' @param color_tax : A taxonomic level to color the points
#' @param taxa : The taxonomic level choose for x-positioning
#' data("GlobalPatterns")
#' res = mt(GlobalPatterns, "SampleType", test="f")
#' plot_mt(res, color_tax = "Phylum") + scale_color_hue()
#' @author Adrien Taudière
#'
#' @return a \code{\link{ggplot}}2 plot of result of a mt test
#'
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
      factor(d$tax, levels = unique(factor(as.character(d[, "Species"]))[rev(order(d$teststat))]))

    p <-
      ggplot(d, aes(x = tax, y = teststat, color = tax_col)) + geom_point(size = 6) +
      theme(axis.text.x = element_text(
        angle = -90,
        hjust = 0,
        vjust = 0.5
      ))
    p
  }
####################################################################################


################################################################################
#' Plot accumulation curves for \code{\link{phyloseq-class}} object
#'
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (required): Name of the factor in physeq@sam_data used to plot different lines
#' @param nbSeq (logical): Either plot accumulation curves using sequences or using samples
#' @param step (integer): distance among points calculated to plot lines. A
#'  low value give better plot but is more time consuming. Only use if nbSeq = TRUE
#' @param by.fact (logical): First merge the OTU table by factor to plot only one line by factor
#' @param ci.col : Color vector for confidence intervall. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to cange plot.
#' @param col : Color vector for lines. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to change plot.
#' @param lwd  (default = 3): thickness for lines. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to cange plot.
#' @param leg (logical): Plot legend or not. Only use if nbSeq = FALSE. If nbSeq = TRUE, you can use ggplot to cange plot.
#' @param printSamplesNames (logical): Print samples names or not? Only use if nbSeq = TRUE.
#' @param CI (default = 2) : Confidence intervall value used to multiply the standard error to plot confidence intervall
#' @param ... Additional arguments passed on to \code{\link{ggplot}} if nbSeq = TRUE
#' or \code{\link{plot}} if nbSeq = FALSE
#'
#' data("GlobalPatterns")
#' GlobalPatterns_Archaea <- subset_taxa(GlobalPatterns,
#'  GlobalPatterns@tax_table[, 1] == 'Archaea')
#' accu_plot(GlobalPatterns_Archaea, 'SampleType', nbSeq = TRUE, by.fact = TRUE)
#'
#' @return A \code{\link{ggplot}}2 plot representing the richness
#' accumulation plot if nbSeq = TRUE, else, if nbSeq = FALSE
#' return a base plot.
#'
#'
#' @author Adrien Taudière
#' @seealso \code{\link[vegan]{specaccum}}
accu_plot <-
  function(physeq,
           fact = NULL,
           nbSeq = TRUE,
           step = NULL,
           by.fact = FALSE,
           ci.col = NULL,
           col = NULL,
           lwd = 3,
           leg = TRUE,
           printSamplesNames = FALSE,
           CI = 2,
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be a phyloseq object")
    }

    if (!nbSeq) {
      factor.interm <-
        eval(parse(text = paste("physeq@sam_data$", fact, sep = "")))
      factor.interm <- as.factor(factor.interm)

      physeq_accu <- as.matrix(t(physeq@otu_table))
      physeq_accu[physeq_accu > 0] <- 1
      accu_all <- vegan::specaccum(physeq_accu)

      accu <- list()
      for (i in 1:nlevels(factor.interm)) {
        accu[[i]] <-
          vegan::specaccum(physeq_accu[factor.interm == levels(factor.interm)[i],])
        #print(paste(round(i/nlevels(factor.interm) * 100), "%"))
      }

      funky.color <-
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
        col <- funky.color(nlevels(factor.interm) + 1)
      }
      if (is.null(ci.col)) {
        transp <- function (col, alpha = 0.5) {
          res <-
            apply(col2rgb(col), 2, function(c)
              rgb(c[1] / 255, c[2] / 255, c[3] / 255, alpha))
          return(res)
        }
        ci.col <- transp(funky.color(nlevels(factor.interm) + 1), 0.3)
      }

      plot(
        accu_all,
        ci.type = "poly",
        ci.col = ci.col[1],
        col = col[1],
        lwd = lwd,
        ci.lty = 0,
        xlab = "Sample",
        ...
      )

      for (i in 1:nlevels(factor.interm)) {
        lines(
          accu[[i]],
          ci.type = "poly",
          ci.col = ci.col[i + 1],
          col = col[i + 1],
          lwd = lwd,
          ci.lty = 0
        )
      }
      if (leg) {
        legend(
          "bottomright",
          c("all", levels(factor.interm)),
          col = col,
          lty = 1,
          lwd = 3
        )
      }
    }

    if (nbSeq) {
      FACT <- as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))

      if (!by.fact) {
        x <- t(physeq@otu_table)
      } else {
        x <- apply(physeq@otu_table, 1, function(x)
          tapply(x, FACT, sum))
      }

      tot <- rowSums(x)
      nr <- nrow(x)

      if (is.null(step)) {
        step = round(max(tot) / 30, 0)
      }
      if (is.null(step)) {
        step = 1
      }

      res <- list()
      n_max <- seq(1, max(tot), by = step)
      out <- lapply(seq_len(nr), function(i) {
        n <- seq(1, tot[i], by = step)
        if (n[length(n)] != tot[i]) {
          n <- c(n, tot[i])
        }
        res_interm <- rarefy(x[i,], n, se = TRUE)
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
          as.factor(unlist(unclass(physeq@sam_data[match(df$.id, sample_names(physeq)), fact])[fact]))
      } else {
        df$fact <- df$.id
      }

      df$ymin <- df$X1 - df$X2 * CI
      df$ymin[is.na(df$ymin)] <- df$X1[is.na(df$ymin)]
      df$ymax <- df$X1 + df$X2 * CI
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
        geom_line() + xlab("Number of sequences") + ylab("Number of OTUs (with standard error)")

      if (printSamplesNames) {
        p + geom_text(data = dff, aes(
          x = xlab ,
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
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (required): Name of the factor to cluster samples by modalities. Need to be in \code{physeq@sam_data}.
#' @param taxa (Default:'Order'): Name of the taxonomic rank of interest
#' @param nbSeq (Default: TRUE): Represent the number of sequences or the number of OTUs (nbSeq = FALSE)
#' @param rarefy (logical): Does each samples modalities need to be rarefy in order to compare them with the same amount of sequences?
#' @param min.prop.tax (Default: 0.01): The minimum proportion for taxon to be ploted
#' @param min.prop.mod (Default: 0.1) : The minimum proportion for modalities to be ploted
#' @param gap.degree : Gap between two neighbour sectors. It can be a single value or a vector. If it is a vector, the first value corresponds to the gap after the first sector.
#' @param start.degree : The starting degree from which the circle begins to draw. Note this degree is measured in the standard polar coordinate which means it is always reverse-clockwise.
#' @param row.col : Color vector for row
#' @param grid.col : Grid colors which correspond to sectors. The length of the vector should be either 1 or the number of sectors. It's preferred that grid.col is a named vector of which names correspond to sectors. If it is not a named vector, the order of grid.col corresponds to order of sectors.
#' @param log10trans (logical): Should sequence be log10 transformed (more precisely by log10(1+x))?
#' @param ... Additional arguments passed on to \code{\link[circlize]{chordDiagram}} or \code{\link[circlize]{circos.par}}
#'
#' data("GlobalPatterns")
#' # GlobalPatterns_Archaea <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == 'Archaea')
#' # otu_circle(GlobalPatterns_Archaea, 'SampleType')
#' @author Adrien Taudière
#'
#' @return A \code{\link[circlize]{chordDiagram}} plot representing the distribution
#' of OTUs or sequences in the different modalities of the factor fact
#'
#' @seealso \code{\link[circlize]{chordDiagram}}
#' @seealso \code{\link[circlize]{circos.par}}

otu_circle <-
  function(physeq = NULL,
           fact = NULL,
           taxa = "Order",
           nbSeq = TRUE,
           rarefy = FALSE,
           min.prop.tax = 0.01,
           min.prop.mod = 0.1,
           gap.degree = NULL,
           start.degree = NULL,
           row.col = NULL,
           grid.col = NULL,
           log10trans = F,
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }

    if (!nbSeq) {
      physeq@otu_table[physeq@otu_table > 0] <- 1
    }

    taxcol <- match(taxa, colnames(physeq@tax_table))
    if (is.na(taxcol)) {
      stop("The taxa argument do not match any taxa rank in physeq@tax_table")
    }

    taxsamp <- match(fact, colnames(physeq@sam_data))
    if (is.na(taxsamp)) {
      stop("The samples argument do not match any sample attributes in physeq@sam_data")
    }

    otu_table_tax <-
      apply(physeq@otu_table, 2, function(x)
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
      print(
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

    # Keep only taxa and modalities with a sufficient proportion (min.prop.tax,
    # min.prop.mod) to plot
    o_t_e_interm <-
      otu_table_ech[(rowSums(otu_table_ech) / sum(otu_table_ech)) > min.prop.mod,
                    (colSums(otu_table_ech) / sum(otu_table_ech)) > min.prop.tax]
    if (nrow(o_t_e_interm) != nrow(otu_table_ech)) {
      print(
        paste(
          "Only ",
          nrow(o_t_e_interm),
          " modalities are plot (",
          round(100 *
                  nrow(o_t_e_interm) /
                  nrow(otu_table_ech), 2),
          "%). Use 'min.prop.mod' to plot more samples.",
          sep = ""
        )
      )
    }

    if (ncol(o_t_e_interm) != ncol(otu_table_ech)) {
      print(
        paste(
          "Only ",
          ncol(o_t_e_interm),
          " taxa are plot (",
          round(100 * ncol(o_t_e_interm) / ncol(otu_table_ech),
                2),
          "%). Use 'min.prop.tax' to plot more taxa",
          sep = ""
        )
      )
    }
    otu_table_ech <- o_t_e_interm

    if (log10trans) {
      otu_table_ech <- apply(otu_table_ech, 2, function(x)
        log10(1 + x))
    }


    if (is.null(gap.degree)) {
      col2keep <- rep(1, ncol(otu_table_ech) - 1)
      row2keep <- rep(1, nrow(otu_table_ech) - 1)
      gap.degree <- c(row2keep, 10, col2keep, 10)

    }
    if (is.null(start.degree)) {
      start.degree <- 170
    }

    funky.color <-
      colorRampPalette(
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

    if (is.null(grid.col)) {
      grid.col <-
        c(funky.color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
    }

    if (is.null(row.col)) {
      row.col <-
        c(funky.color(nrow(otu_table_ech)), rep("grey", ncol(otu_table_ech)))
    }

    circlize::circos.par(gap.degree = gap.degree, start.degree = start.degree, ...)
    circlize::chordDiagram(otu_table_ech, row.col = row.col, grid.col = grid.col, ...)
    circlize::circos.clear()
  }
################################################################################

################################################################################
#' Sankey plot of \code{\link{phyloseq-class}} object
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (Optional): Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}.
#' @param taxa (Default: c(1:4)): a vector of taxonomic rank to plot
#' @param nbSeq (Default: FALSE): Represent the number of sequences or the number of OTUs (nbSeq = FALSE). Note that ploting the number of sequences is slower.
#' @param min.prop.tax (Default: 0): The minimum proportion for taxon to be ploted. EXPERIMENTAL. For the moment each links below the min.prop. tax is discard from the sankey network resulting in sometimes weird plot.
#' @param tax2remove : a vector of taxonomic groups to remove from the analysis (e.g. \code{c('Incertae sedis', 'unidentified')})
#' @param units : character string describing physical units (if any) for Value
#' @param Symbol2sub (default = c('\\.', '-')): vector of symbol to delete in the taxonomy
#' @param ... Additional arguments passed on to \code{\link[networkD3]{sankeyNetwork}}
#'
#' data("GlobalPatterns")
#' #GlobalPatterns_Archaea <- subset_taxa(GlobalPatterns, GlobalPatterns@tax_table[, 1] == 'Archaea')
#' #sankey_phyloseq(GlobalPatterns_Archaea)
#' #sankey_phyloseq(GlobalPatterns_Archaea, fact = 'SampleType')
#' sankey_phyloseq(GlobalPatterns, taxa = c(1:5), min.prop.tax = 0.01)
#' sankey_phyloseq(GlobalPatterns, taxa = c(2:6), min.prop.tax = 0.01, nbSeq = TRUE)
#' @author Adrien Taudière
#'
#' @return A \code{\link[networkD3]{sankeyNetwork}} plot representing the taxonomic distribution
#' of OTUs or sequences. If \code{fact} is set, represent the distribution of
#' the last taxonomic level in the modalities of \code{fact}
#'
#' @seealso \code{\link[networkD3]{sankeyNetwork}}

sankey_phyloseq <-
  function(physeq = NULL,
           fact = NULL,
           taxa = c(1:4),
           nbSeq = FALSE,
           min.prop.tax = 0,
           tax2remove = NULL,
           units = NULL,
           Symbol2sub = c("\\.", "-"),
           ...) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }

    if (!nbSeq) {
      physeq@otu_table[physeq@otu_table > 0] <- 1
      mat.interm <- matrix()
      mat <- matrix(ncol = 3)
      colnames (mat) <- c("Var1", "Var2", "value")
      for (i in 1:(length(taxa) - 1)) {
        res.interm <-
          table(physeq@tax_table[, taxa[i]], physeq@tax_table[, taxa[i + 1]])
        mat.interm <- reshape2::melt(res.interm)
        mat.interm <- mat.interm[mat.interm[, 3] > 0, ]
        mat <- rbind(mat, mat.interm)
      }
    } else if (nbSeq) {
      mat.interm <- matrix()
      mat <- matrix(ncol = 3)
      colnames (mat) <- c("Var1", "Var2", "value")
      tax_table.interm <-
        physeq@tax_table[rep(1:dim(physeq@tax_table)[1], times = taxa_sums(physeq))]

      for (i in 1:(length(taxa) - 1)) {
        res.interm <-
          table(tax_table.interm[, taxa[i]], tax_table.interm[, taxa[i + 1]])
        mat.interm <- reshape2::melt(res.interm)
        mat.interm <- mat.interm[mat.interm[, 3] > 0, ]
        mat <- rbind(mat, mat.interm)
      }
    }

    if (!is.null(fact)) {
      NetMatrix2Links <- function(m = NULL) {
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

      mat.interm <-
        apply(physeq@otu_table, 1, function(x)
          tapply(x, physeq@sam_data[, fact],
                 sum))

      if (!nbSeq) {
        mat.interm <-
          apply(mat.interm, 1, function(x)
            tapply(x, physeq@tax_table[,
                                       taxa[length(taxa)]], function(x)
                                         sum(x > 0)))
      } else if (nbSeq) {
        mat.interm <-
          apply(mat.interm, 1, function(x)
            tapply(x, physeq@tax_table[,
                                       taxa[length(taxa)]], sum))
      }

      sampLinks <- NetMatrix2Links(mat.interm)
      sampLinks[, 2] <- toupper(sampLinks[, 2])
      colnames(sampLinks) <- colnames(mat)
      mat <- rbind(mat, sampLinks)
    }

    mat <- as.data.frame(mat[rowSums(is.na(mat)) == 0,])
    mat[, 3] <- as.numeric(as.vector(mat[, 3]))
    mat <- mat[rowSums(is.na(mat)) == 0,]


    if (!is.null(tax2remove)) {
      mat <- mat[!mat[, 1] %in% tax2remove,]
      mat <- mat[!mat[, 2] %in% tax2remove,]
    }

    if (min.prop.tax != 0) {
      min.nb.tax <- min.prop.tax * sum(mat[, 3]) / length(taxa)
      mat <- mat[mat[, 3] >= min.nb.tax,]
    }

    for (i in 1:length(Symbol2sub)) {
      mat <- apply(mat, 2, function(x)
        gsub(Symbol2sub[i], "", x))
    }

    taxSank <- list()
    namesNodes <- unique(c(as.vector(mat[, 1]), as.vector(mat[, 2])))
    namesNodes <- namesNodes[!is.na(namesNodes)]
    taxSank$nodes <-
      data.frame((1:length(namesNodes)) - 1, namesNodes)
    names(taxSank$nodes) <- c("code", "name")
    mat2 <- mat
    for (i in 1:nrow(taxSank$nodes)) {
      mat2[, 1] <-
        gsub(paste("\\<", taxSank$nodes[i, 2], "\\>", sep = ""),
             taxSank$nodes[i,
                           1],
             mat2[, 1])
      mat2[, 2] <-
        gsub(paste("\\<", taxSank$nodes[i, 2], "\\>", sep = ""),
             taxSank$nodes[i,
                           1],
             mat2[, 2])
    }

    taxSank$links <- apply(mat2, 2, as.numeric)
    taxSank$links <-
      data.frame(taxSank$links[rowSums(is.na(taxSank$links)) == 0,])
    taxSank$nodes <- as.data.frame(as.character(taxSank$nodes[, 2]))
    names(taxSank$nodes) <- c("name")
    names(taxSank$links) <- c("source", "target", "value")
    if (is.null(units)) {
      if (!nbSeq) {
        units <- "OTUs"
      } else if (nbSeq) {
        units <- "Sequences"
      }
    }
    networkD3::sankeyNetwork(
      Links = taxSank$links,
      Nodes = taxSank$nodes,
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
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (required): Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}.
#' @param min.nb.seq (Default: 0)): minimum number of sequences by OTUs by samples
#' to take into count this OTUs in this sample. For example, if min.nb.seq=2, each value of 2 or less in the OTU table
#' will be change into 0 for the analysis
#' @param printValues (logical) : Print (or not) the table of number of OTUs
#' for each combination. If printValues is TRUE the object is not a ggplot object.
#' Please use printValues = FALSE if you want to add ggplot function (cf example).
#'
#'
#' data("enterotype")
#' venn_phyloseq(enterotype, fact = 'SeqTech')
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus')
#' venn_phyloseq(enterotype, fact = 'Nationality', printValues = F)
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus', printValues = F) + scale_fill_hue()
#' venn_phyloseq(enterotype, fact = 'ClinicalStatus', printValues = F) + scale_fill_hue()
#'
#' @return A \code{\link{ggplot}}2 plot representing Venn diagramm of
#' modalities of the argument \code{factor}
#'
#' @author Adrien Taudière
#' @seealso \code{\link[venneuler]{venneuler}}

venn_phyloseq <-
  function(physeq,
           fact,
           min.nb.seq = 0,
           printValues = TRUE) {
    if (!inherits(physeq, "phyloseq")) {
      stop("physeq must be an object of class 'phyloseq'")
    }

    moda <- as.factor(unlist(unclass(physeq@sam_data[, fact])[fact]))
    if (length(moda) != dim(physeq@otu_table)[1]) {
      data_Venn <-
        t(apply(physeq@otu_table, 1, function(x)
          by(x, moda, max)))
    } else if (length(moda) != dim(physeq@otu_table)[2]) {
      data_Venn <-
        t(apply(t(physeq@otu_table), 1, function(x)
          by(x, moda, max)))
    } else {
      stop("The factor length and the number of samples must be identical")
    }
    combinations <- data_Venn > min.nb.seq

    e <- new.env(TRUE, emptyenv())
    cn <- colnames(combinations)
    for (i in seq.int(dim(combinations)[1]))
      if (any(combinations[i,])) {
        ec <- paste(cn[combinations[i,]], collapse = "&")
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

    VENN <- venneuler::venneuler(data_Venn > min.nb.seq)
    venn_res <-
      data.frame(
        x = VENN$centers[, 1],
        y = VENN$centers[, 2],
        radius = VENN$diameters / 2
      )

    nmod <- nrow(venn_res)
    x1 <- list()
    for (i in 1:nmod) {
      x1[[i]] <- grep(rownames(venn_res)[i], table_value$combinations)
    }

    for (i in 1:nrow(table_value)) {
      table_value$x[i] <- mean(VENN$centers[, "x"][unlist(lapply(x1,
                                                                 function(x)
                                                                   sum(x %in% i) > 0))])
      table_value$y[i] <- mean(VENN$centers[, "y"][unlist(lapply(x1,
                                                                 function(x)
                                                                   sum(x %in% i) > 0))])
    }

    df <- venn_res
    df$xlab <- df$x + (df$x - mean(df$x))
    df$ylab <- df$y + (df$y - mean(df$y))

    circularise <- function(d, n = 360) {
      angle <- seq(-pi, pi, length = n)
      make_circle <- function(x, y, r, Modality) {
        data.frame(x = x + r * cos(angle),
                   y = y + r * sin(angle),
                   Modality)
      }
      lmat <- mapply(
        make_circle,
        Modality = rownames(d),
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
                              aes(x, y, group = Modality, fill = Modality),
                              alpha = 0.5) + theme_void()

    if (printValues) {
      g_legend <- function(a.gplot) {
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
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
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <-
  function(...,
           plotlist = NULL,
           file,
           cols = 1,
           layout = NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
                       ncol = cols,
                       nrow = ceiling(numPlots / cols))
    }

    if (numPlots == 1) {
      print(plots[[1]])

    } else {
      # Set up the page
      grid::grid.newpage()
      grid::pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

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
#'
#' @param physeq (required): A \code{\link{phyloseq-class}} object
#' @param variable (required): The variable to test
#' @param color_fac (optional): The variable to color the barplot
#'
#' @return A list of 4 ggplot2 plot.
#'   - plot_Hill_0 : the boxplot of Hill number 0 (= species richness) against the variable
#'   - plot_Hill_1 : the boxplot of Hill number 1 (= Shannon index) against the variable
#'   - plot_Hill_2 : the boxplot of Hill number 2 (= Simpson index) against the variable
#'   - plot_tuckey : plot the result of the Tuckey HSD test
#'
#' @export
#'
hill_phyloseq <- function(physeq, variable, color_fac) {
  var <- sym(variable)
  color_fac <- sym(color_fac)

  otuHill <-
    vegan::renyi(physeq@otu_table, scale = c(0, 1, 2), hill = T)
  colnames(otuHill) <- c("Hill_0", "Hill_1", "Hill_2")

  df_hill <- data.frame(otuHill, physeq@sam_data)
  df_hill[, c(1:3)] <- apply(df_hill[, c(1:3)], 2, as.numeric)

  p_0 <- ggplot(df_hill, aes(group = !!var, Hill_0))  +
    geom_boxplot(outlier.size = 2, aes(colour = as.factor(!!color_fac)))
  p_1 <- ggplot(df_hill, aes(group = !!var, Hill_1))  +
    geom_boxplot(outlier.size = 2, aes(colour = as.factor(!!color_fac)))
  p_2 <- ggplot(df_hill, aes(group = !!var, Hill_2))  +
    geom_boxplot(outlier.size = 2, aes(colour = as.factor(!!color_fac)))

  p_var <- hill_tuckey_phyloseq(physeq, variable)

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
#'
#' @param physeq (required): A \code{\link{phyloseq-class}} object
#' @param abundance (Default: TRUE): Does the number of sequences is print
#' @param ... Other argument for the datatable function
#'
#' @author Adrien Taudière
#' @return a datatable
#' @export
#'
#' @example
#' data("GlobalPatterns")
#' tax_datatable(subset_taxa(GlobalPatterns, rowSums(GlobalPatterns@otu_table)>10000))
tax_datatable <- function(physeq, abundance = TRUE, ...) {
  df <- as.data.frame(physeq@tax_table)

  if (abundance) {
    if (physeq@otu_table@taxa_are_rows) {
      df$nb_seq <- rowSums(physeq@otu_table)
    } else {
      df$nb_seq <- colSums(physeq@otu_table)
    }
  }

  DT::datatable(df, ...) %>%  DT::formatStyle(
    'nb_seq',
    background = DT::styleColorBar(df$nb_seq, 'steelblue'),
    backgroundSize = '100% 90%',
    backgroundRepeat = 'no-repeat',
    backgroundPosition = 'center'
  ) %>%
    DT::formatStyle(
      'Order',
      borderBottom = '5px solid',
      # textDecorationLine = 'underline',
      borderColor = DT::styleEqual(unique(df$Order), viridis::viridis(length(unique(
        df$Order
      ))))
    )
}
################################################################################

################################################################################
#' Summarise a \code{\link{phyloseq-class}} object using a plot.
#'
#' @param physeq (required): A \code{\link{phyloseq-class}} object
#'
#' @return A ggplot2 object
#' @export
summary_plot_phyloseq <- function(physeq) {
  d = data.frame(
    y1 = c(1, 1, 5.3, 1),
    y2 = c(5, 5, 7.5, 5),
    x1 = c(1, 3.15, 1, 4.3),
    x2 = c(3, 4.15, 3, 5.5),
    names = c('@otu_table', '@tax_table', '@sam_data', '@refseq'),
    nb_values = c(
      paste(
        format(ncol(physeq@otu_table), big.mark = " "),
        "taxa\n",
        format(sum(physeq@otu_table), big.mark = " "),
        "sequences\n",
        format(sum(physeq@otu_table > 0), big.mark = " "),
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
              nchar(physeq@refseq)
            ), 2), "+/-",
            round(sd(
              nchar(physeq@refseq)
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
