################################################################################
#' @title Performs graph-based permutation tests on phyloseq object
#' @description
#' `r lifecycle::badge("maturing")`
#'
#' A wrapper of [phyloseqGraphTest::graph_perm_test()] for quick plot with
#' important statistics

#' @inheritParams clean_pq
#' @param fact (required) Name of the factor to cluster samples by modalities.
#'   Need to be in \code{physeq@sam_data}. This should be a factor
#'   with two or more levels.
#' @param merge_sample_by a vector to determine
#'   which samples to merge using [merge_samples2()] function.
#'   Need to be in \code{physeq@sam_data}
#' @param nperm (int) The number of permutations to perform.
#' @param return_plot (logical) Do we return only the result
#'   of the test or do we plot the result?
#' @param title The title of the Graph.
#' @param na_remove (logical, default FALSE) If set to TRUE, remove samples with
#'   NA in the variables set in formula.
#' @param ... Other params for be passed on to
#'   [phyloseqGraphTest::graph_perm_test()] function
#'
#' @examples
#' \donttest{
#' if (requireNamespace("phyloseqGraphTest")) {
#'   data(enterotype)
#'   graph_test_pq(enterotype, fact = "SeqTech")
#'   graph_test_pq(enterotype, fact = "Enterotype", na_remove = TRUE)
#' }
#' }
#' @author Adrien Taudière
#'
#' @return A \code{\link{ggplot}}2 plot with a subtitle indicating the pvalue
#' and the number of permutations
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please cite `phyloseqGraphTest` package.
#' @export

graph_test_pq <- function(physeq,
                          fact,
                          merge_sample_by = NULL,
                          nperm = 999,
                          return_plot = TRUE,
                          title = "Graph Test",
                          na_remove = FALSE,
                          ...) {
  verify_pq(physeq)

  if (!is.null(merge_sample_by)) {
    physeq <- merge_samples2(physeq, merge_sample_by)
    physeq <- clean_pq(physeq)
  }

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

  res_graph_test <- phyloseqGraphTest::graph_perm_test(physeq,
    sampletype = fact,
    nperm = nperm,
    ...
  )
  if (!return_plot) {
    return(res_graph_test)
  } else {
    p <- phyloseqGraphTest::plot_test_network(res_graph_test) +
      labs(
        title = title,
        subtitle = paste(
          "pvalue = ",
          res_graph_test$pval,
          "(",
          length(res_graph_test$perm),
          " permutations)"
        )
      )
    return(p)
  }
}
################################################################################

################################################################################
#' @title Permanova on a phyloseq object
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [vegan::adonis2()] function in the case of `physeq` object.
#' @inheritParams clean_pq
#' @param formula (required) the right part of a formula for [vegan::adonis2()].
#'   Variables must be present in the `physeq@sam_data` slot.
#' @param dist_method (default "bray") the distance used. See
#'   [phyloseq::distance()] for all available distances or run
#'   [phyloseq::distanceMethodList()].
#'   For aitchison and robust.aitchison distance, [vegan::vegdist()]
#'   function is directly used.
#' @param merge_sample_by a vector to determine
#'   which samples to merge using the [merge_samples2()]
#'   function. Need to be in `physeq@sam_data`
#' @param na_remove (logical, default FALSE) If set to TRUE, remove samples with
#'   NA in the variables set in formula.
#' @param correction_for_sample_size (logical, default FALSE) If set to TRUE,
#'   the sample size (number of sequences by samples) is add to formula in
#'   the form `y~Library_Size + Biological_Effect` following recommendation of
#'   [Weiss et al. 2017](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y).
#'   `correction_for_sample_size` overcome `rarefy_nb_seqs` if both are TRUE.
#' @param rarefy_nb_seqs (logical, default FALSE) Rarefy each sample
#'   (before merging if merge_sample_by is set) using
#'    `phyloseq::rarefy_even_depth()`.
#'   if `correction_for_sample_size` is TRUE, rarefy_nb_seqs will have no
#'   effect.
#' @param verbose (logical, default TRUE) If TRUE, prompt some messages.
#' @param ... Other arguments passed on to [vegan::adonis2()] function.
#' @return The function returns an anova.cca result object with a
#'   new column for partial R^2. See help of [vegan::adonis2()] for
#'   more information.
#' @examples
#' data(enterotype)
#' \donttest{
#' adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE)
#' adonis_pq(enterotype, "SeqTech", dist_method = "jaccard")
#' adonis_pq(enterotype, "SeqTech", dist_method = "robust.aitchison")
#' }
#' @export
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `vegan::adonis2()` if you
#'   use this function.

adonis_pq <- function(physeq,
                      formula,
                      dist_method = "bray",
                      merge_sample_by = NULL,
                      na_remove = FALSE,
                      correction_for_sample_size = FALSE,
                      rarefy_nb_seqs = FALSE,
                      verbose = TRUE,
                      ...) {
  physeq <- taxa_as_columns(physeq)

  if (dist_method %in% c("aitchison", "robust.aitchison")) {
    phy_dist <-
      paste0(
        'vegan::vegdist(as.matrix(physeq@otu_table), method="',
        dist_method,
        '")'
      )
  } else {
    phy_dist <-
      paste0('phyloseq::distance(physeq, method="', dist_method, '")')
  }

  .formula <- stats::reformulate(formula, response = phy_dist)
  termf <- stats::terms(.formula)
  term_lab <- attr(termf, "term.labels")[attr(termf, "order") == 1]

  verify_pq(physeq)

  if (na_remove) {
    new_physeq <- physeq
    for (tl in term_lab) {
      if (verbose) {
        message(tl)
      }
      new_physeq <-
        subset_samples_pq(new_physeq, !is.na(physeq@sam_data[[tl]]))
    }
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

  if (!is.null(merge_sample_by)) {
    physeq <- merge_samples2(physeq, merge_sample_by)
    physeq <- clean_pq(physeq)
  }

  if (correction_for_sample_size) {
    formula <- paste0("sample_size+", formula)
    .formula <- stats::reformulate(formula, response = phy_dist)
  } else if (rarefy_nb_seqs) {
    physeq <- rarefy_even_depth(physeq)
    physeq <- clean_pq(physeq)
  }

  verify_pq(physeq)
  metadata <- as(sample_data(physeq), "data.frame")
  if (correction_for_sample_size) {
    metadata$sample_size <- sample_sums(physeq)
  }

  res_ado <- vegan::adonis2(.formula, data = metadata, ...)
  return(res_ado)
}
################################################################################

################################################################################
#' Permanova (adonis) on permutations of rarefaction even depth
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams adonis_pq
#' @param nperm (int) The number of permutations to perform.
#' @param progress_bar (logical, default TRUE) Do we print progress during
#'   the calculation.
#' @param quantile_prob (float, `[0:1]`) the value to compute the quantile.
#'   Minimum quantile is compute using 1-quantile_prob.
#' @param sample.size (int) A single integer value equal to the number of
#'   reads being simulated, also known as the depth. See
#'   [phyloseq::rarefy_even_depth()].
#' @param ... Other params for be passed on to [adonis_pq()] function
#'
#' @return A list of three dataframe representing the mean, the minimum quantile
#'  and the maximum quantile value for adonis results. See [adonis_pq()].
#' @export
#' @author Adrien Taudière
#' @seealso [adonis_pq()]
#' @examples
#' if (requireNamespace("vegan")) {
#'   data_fungi_woNA <-
#'     subset_samples(data_fungi, !is.na(Time) & !is.na(Height))
#'   adonis_rarperm_pq(data_fungi_woNA, "Time*Height", na_remove = TRUE, nperm = 3)
#' }
adonis_rarperm_pq <- function(physeq,
                              formula,
                              dist_method = "bray",
                              merge_sample_by = NULL,
                              na_remove = FALSE,
                              rarefy_nb_seqs = FALSE,
                              verbose = TRUE,
                              nperm = 99,
                              progress_bar = TRUE,
                              quantile_prob = 0.975,
                              sample.size = min(sample_sums(physeq)),
                              ...) {
  res_perm <- list()
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
    res_perm[[i]] <-
      adonis_pq(
        rarefy_even_depth(
          physeq,
          rngseed = i,
          sample.size = sample.size,
          verbose = verbose
        ),
        formula,
        dist_method = dist_method,
        merge_sample_by = merge_sample_by,
        na_remove = na_remove,
        correction_for_sample_size = FALSE,
        rarefy_nb_seqs = rarefy_nb_seqs,
        verbose = verbose,
        sample.size = sample.size,
        ...
      )
    if (progress_bar) {
      setTxtProgressBar(pb, i)
    }
  }
  res_adonis <- list()
  res_adonis[["mean"]] <-
    apply(array(unlist(res_perm), c(dim(
      as.data.frame(res_perm[[1]])
    ), nperm)), c(1, 2), mean)
  colnames(res_adonis[["mean"]]) <- colnames(res_perm[[1]])
  rownames(res_adonis[["mean"]]) <- rownames(res_perm[[1]])


  res_adonis[["quantile_min"]] <-
    apply(
      array(unlist(res_perm), c(dim(
        as.data.frame(res_perm[[1]])
      ), nperm)),
      c(1, 2),
      quantile,
      na.rm = TRUE,
      probs = 1 - quantile_prob
    )
  colnames(res_adonis[["quantile_min"]]) <- colnames(res_perm[[1]])
  rownames(res_adonis[["quantile_min"]]) <- rownames(res_perm[[1]])

  res_adonis[["quantile_max"]] <-
    apply(
      array(unlist(res_perm), c(dim(
        as.data.frame(res_perm[[1]])
      ), nperm)),
      c(1, 2),
      quantile,
      na.rm = TRUE,
      probs = quantile_prob
    )
  colnames(res_adonis[["quantile_max"]]) <- colnames(res_perm[[1]])
  rownames(res_adonis[["quantile_max"]]) <- rownames(res_perm[[1]])

  return(res_adonis)
}
################################################################################

################################################################################
#' @title Compute and test local contributions to beta diversity (LCBD) of
#'   samples
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [adespatial::beta.div()] function in the case of `physeq`
#'   object.
#' @inheritParams clean_pq
#'
#' @param p_adjust_method (chr, default "BH"): the method used to adjust p-value
#' @param ... Others arguments passed on to [adespatial::beta.div()] function
#'
#' @return An object of class `beta.div` see [adespatial::beta.div()] function
#'   for more information
#' @export
#' @seealso [plot_LCBD_pq], [adespatial::beta.div()]
#' @examples
#' if (requireNamespace("adespatial")) {
#'   res <- LCBD_pq(data_fungi_sp_known, nperm = 5)
#'   str(res)
#'   length(res$LCBD)
#'   length(res$SCBD)
#' }
#' \donttest{
#' if (requireNamespace("adespatial")) {
#'   LCBD_pq(data_fungi_sp_known, nperm = 5, method = "jaccard")
#' }
#' }
#'
#' @author Adrien Taudière
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `adespatial::beta.div()` if you
#'   use this function.
LCBD_pq <- function(physeq,
                    p_adjust_method = "BH",
                    ...) {
  physeq <- taxa_as_columns(physeq)

  mat <- as.matrix(unclass(physeq@otu_table))
  resBeta <- adespatial::beta.div(mat, adj = FALSE, ...)

  resBeta$p.adj <-
    p.adjust(resBeta$p.LCBD, method = p_adjust_method)
  return(resBeta)
}
################################################################################

################################################################################
#' @title Plot and test local contributions to beta diversity (LCBD) of samples
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [adespatial::beta.div()] function in the case of `physeq`
#'   object.
#' @inheritParams clean_pq
#'
#' @param p_adjust_method (chr, default "BH"): the method used to adjust p-value
#' @param pval (int, default 0.05): the value to determine the significance of
#'   LCBD
#' @param sam_variables A vector of variables names present in the `sam_data`
#'   slot to plot alongside the LCBD value
#' @param only_plot_significant (logical, default TRUE) Do we plot all LCBD
#'   values or only the significant ones
#' @param ... Others arguments passed on to [adespatial::beta.div()] function
#'
#' @return A ggplot2 object build with the package patchwork
#' @export
#' @seealso [LCBD_pq], [adespatial::beta.div()]
#'
#' @examples
#' data(data_fungi)
#' if (requireNamespace("adespatial")) {
#'   plot_LCBD_pq(data_fungi_mini,
#'     nperm = 100, only_plot_significant = FALSE,
#'     pval = 0.2
#'   )
#' }
#' \donttest{
#' if (requireNamespace("adespatial")) {
#'   plot_LCBD_pq(data_fungi_mini,
#'     nperm = 100, only_plot_significant = TRUE,
#'     pval = 0.2
#'   )
#'   if (requireNamespace("patchwork")) {
#'     plot_LCBD_pq(data_fungi_mini,
#'       nperm = 100, only_plot_significant = FALSE,
#'       sam_variables = c("Time", "Height")
#'     )
#'     plot_LCBD_pq(data_fungi_mini,
#'       nperm = 100, only_plot_significant = TRUE, pval = 0.2,
#'       sam_variables = c("Time", "Height", "Tree_name")
#'     ) &
#'       theme(
#'         legend.key.size = unit(0.4, "cm"),
#'         legend.text = element_text(size = 10),
#'         axis.title.x = element_text(size = 6)
#'       )
#'   }
#' }
#' }
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `vegan::beta.div()` if you
#'   use this function.
plot_LCBD_pq <- function(physeq,
                         p_adjust_method = "BH",
                         pval = 0.05,
                         sam_variables = NULL,
                         only_plot_significant = TRUE,
                         ...) {
  resBeta <- LCBD_pq(physeq,
    p_adjust_method = p_adjust_method,
    ...
  )

  sam_data <- data.frame(physeq@sam_data)

  if (sum(is.na(match(
    rownames(sam_data), names(resBeta$LCBD)
  ))) > 0) {
    warning("At least one sample was removed by the beta.div function")
    sam_data <-
      sam_data[rownames(sam_data) %in% names(resBeta$LCBD), ]
  }

  resLCBD <- tibble(
    "LCBD" = resBeta$LCBD,
    "p.LCBD" = resBeta$p.LCBD,
    "p.adj" = resBeta$p.adj,
    sam_data
  )

  if (only_plot_significant) {
    p_LCBD <-
      ggplot(
        filter(resLCBD, p.adj < pval),
        aes(
          x = LCBD,
          y = reorder(Sample_names, -LCBD, sum)
        )
      ) +
      geom_point(size = 4)
    if (is.null(sam_variables)) {
      return(p_LCBD)
    } else {
      p_heatmap <- list()
      for (i in seq_len(length(sam_variables))) {
        p_heatmap[[i]] <- ggplot(filter(resLCBD, p.adj < pval)) +
          geom_tile(inherit.aes = FALSE, aes(
            y = reorder(Sample_names, -LCBD, sum),
            x = i,
            fill = .data[[sam_variables[[i]]]]
          )) +
          theme_void() +
          theme(axis.title.x = element_text()) +
          xlab(sam_variables[[i]])
      }

      p <-
        p_LCBD + patchwork::wrap_plots(p_heatmap) + patchwork::plot_layout(
          widths = c(3, 1), guides =
            "collect"
        )
      return(p)
    }
  } else {
    p_LCBD <-
      ggplot(resLCBD, aes(x = LCBD, y = reorder(Sample_names, -LCBD, sum))) +
      geom_point(aes(color = p.adj < pval))
    if (is.null(sam_variables)) {
      return(p_LCBD)
    } else {
      p_heatmap <- list()
      for (i in seq_len(length(sam_variables))) {
        p_heatmap[[i]] <- ggplot(resLCBD) +
          geom_tile(inherit.aes = FALSE, aes(
            y = reorder(Sample_names, -LCBD, sum),
            x = i,
            fill = .data[[sam_variables[[i]]]]
          )) +
          theme_void() +
          theme(axis.title.x = element_text()) +
          xlab(sam_variables[[i]])
      }

      p <-
        p_LCBD + patchwork::wrap_plots(p_heatmap) + patchwork::plot_layout(
          widths = c(3, 1), guides =
            "collect"
        )
      return(p)
    }
  }
}
################################################################################

################################################################################
#' @title Plot species contributions to beta diversity (SCBD) of samples
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [adespatial::beta.div()] function in the case of `physeq`
#'   object.
#' @inheritParams clean_pq
#' @param tax_level Taxonomic level to used in y axis
#' @param tax_col Taxonomic level to colored points
#' @param min_SCBD (default 0.01) the minimum SCBD value
#'   to plot the taxa

#' @param ... Others arguments passed on to [adespatial::beta.div()] function
#'
#' @return A ggplot2 object build with the package patchwork
#' @export
#' @seealso [LCBD_pq], [adespatial::beta.div()]
#'
#' @examples
#' data(data_fungi)
#' if (requireNamespace("adespatial")) {
#'   plot_SCBD_pq(data_fungi) +
#'     geom_text(aes(label = paste(Genus, Species)), hjust = 1, vjust = 2) +
#'     xlim(c(0, NA))
#' }
#' \donttest{
#' if (requireNamespace("adespatial")) {
#'   plot_SCBD_pq(data_fungi, tax_level = "Class", tax_col = "Phylum", min_SCBD = 0) +
#'     geom_jitter()
#' }
#' }
#'
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `vegan::beta.div()` if you
#'   use this function.
plot_SCBD_pq <- function(physeq,
                         tax_level = "ASV",
                         tax_col = "Order",
                         min_SCBD = 0.01,
                         ...) {
  resBeta <- LCBD_pq(physeq, nperm = 0, ...)

  tax_tab <- data.frame(physeq@tax_table)

  resSCBD <- tibble(
    "ASV" = taxa_names(physeq),
    "SCBD" = resBeta$SCBD,
    tax_tab
  )

  p_SCBD <- ggplot(
    filter(resSCBD, SCBD > min_SCBD),
    aes(
      x = SCBD, y = factor(.data[[tax_level]]),
      color = .data[[tax_col]]
    )
  ) +
    geom_point()

  return(p_SCBD)
}
################################################################################

################################################################################
#' @title Test and plot multipatt result
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [indicspecies::multipatt()] function in the case of
#'   `physeq` object.
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot
#'    different lines
#' @param p_adjust_method (chr, default "BH"): the method used to adjust p-value
#' @param pval (int, default 0.05): the value to determine the significance of
#'   LCBD
#' @param control see `?indicspecies::multipatt()`
#' @param ... Others arguments passed on to [indicspecies::multipatt()] function
#'
#' @return A ggplot2 object
#' @export
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' if (requireNamespace("indicspecies")) {
#'   data(data_fungi)
#'   data_fungi_ab <- subset_taxa_pq(data_fungi, taxa_sums(data_fungi) > 10000)
#'   multipatt_pq(subset_samples(data_fungi_ab, !is.na(Time)), fact = "Time")
#' }
#' \donttest{
#' if (requireNamespace("indicspecies")) {
#'   multipatt_pq(subset_samples(data_fungi_ab, !is.na(Time)),
#'     fact = "Time",
#'     max.order = 1, control = permute::how(nperm = 99)
#'   )
#' }
#' }
#'
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `indicspecies::multipatt()` if you
#'   use this function.

multipatt_pq <- function(physeq,
                         fact,
                         p_adjust_method = "BH",
                         pval = 0.05,
                         control = permute::how(nperm = 999),
                         ...) {
  physeq <- taxa_as_columns(physeq)

  res <-
    indicspecies::multipatt(as.matrix(physeq@otu_table),
      physeq@sam_data[[fact]],
      control = control,
      ...
    )

  res_df <- res$sign
  res_df$p.adj <- p.adjust(res_df$p.value, method = p_adjust_method)
  res_df$ASV_names <- rownames(res_df)
  res_df_signif <-
    res_df %>%
    filter(p.adj < pval) %>%
    tidyr::pivot_longer(cols = starts_with("s."))

  p <- ggplot(
    res_df_signif,
    aes(
      x = ASV_names,
      y = name,
      size = 2 * value,
      color = stat
    )
  ) +
    geom_point() +
    theme(axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ))
  return(p)
}
################################################################################

################################################################################
#' Run ANCOMBC2 on phyloseq object
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [ANCOMBC::ancombc2()] function
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot
#'    different lines
#' @param levels_fact (default NULL) The order of the level in the factor.
#'   Used for reorder levels and select levels (filter out levels not present
#'   en levels_fact)
#' @param tax_level The taxonomic level passed on to [ANCOMBC::ancombc2()]
#' @param ... Other arguments passed on to [ANCOMBC::ancombc2()] function.
#' @return The result of [ANCOMBC::ancombc2()] function
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mia")) {
#'   data_fungi_mini@tax_table <- phyloseq::tax_table(cbind(
#'     data_fungi_mini@tax_table,
#'     "taxon" = taxa_names(data_fungi_mini)
#'   ))
#'   res_height <- ancombc_pq(
#'     data_fungi_mini,
#'     fact = "Height",
#'     levels_fact = c("Low", "High"),
#'     verbose = TRUE
#'   )
#'
#'   ggplot(
#'     res_height$res,
#'     aes(
#'       y = reorder(taxon, lfc_HeightHigh),
#'       x = lfc_HeightHigh,
#'       color = diff_HeightHigh
#'     )
#'   ) +
#'     geom_vline(xintercept = 0) +
#'     geom_segment(aes(
#'       xend = 0, y = reorder(taxon, lfc_HeightHigh),
#'       yend = reorder(taxon, lfc_HeightHigh)
#'     ), color = "darkgrey") +
#'     geom_point()
#'
#'   res_time <- ancombc_pq(
#'     data_fungi_mini,
#'     fact = "Time",
#'     levels_fact = c("0", "15"),
#'     tax_level = "Family",
#'     verbose = TRUE
#'   )
#' }
#' }
#' @author Adrien Taudière
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `ANCOMBC::ancombc2()` if you
#'   use this function.
ancombc_pq <- function(physeq, fact, levels_fact = NULL, tax_level = "Class", ...) {
  if (!is.null(levels_fact)) {
    physeq <- subset_samples_pq(physeq, as.vector(physeq@sam_data[, fact])[[1]] %in% levels_fact)
  }
  tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(physeq)
  if (!is.null(levels_fact)) {
    SummarizedExperiment::colData(tse)[[fact]] <- factor(tse[[fact]], levels = levels_fact)
  }
  res_ancomb <- ANCOMBC::ancombc2(
    data = tse,
    assay_name = "counts",
    fix_formula = fact,
    tax_level = tax_level,
    group = fact,
    ...
  )
  return(res_ancomb)
}
################################################################################


################################################################################
#' Filter ancombc_pq results
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @param ancombc_res (required) the result of the ancombc_pq function
#'   For the moment only bimodal factors are possible.
#' @param filter_passed (logical, default TRUE) Do we filter using the column
#'   passed_ss? The passed_ss value is  TRUE if the taxon passed the sensitivity
#'   analysis, i.e., adding different pseudo-counts to 0s would not change the results.
#' @param filter_diff (logical, default TRUE) Do we filter using the column
#'   diff? The diff value is TRUE if the taxon is significant
#'   (has q less than alpha)
#' @param min_abs_lfc (integer, default0) Minimum absolute value to filter
#'   results based on Log Fold Change. For ex. a value of 1 filter out taxa
#'   for which the abundance in a given level of the modality is not at least
#'   the double of the abundance in the other level.
#'
#' @return A data.frame with the same number of columns than the `ancombc_res`
#'   param but with less (or equal) numbers of rows
#' @export
#'
#' @seealso [ancombc_pq()], [plot_ancombc_pq()]
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mia")) {
#'   data_fungi_mini@tax_table <- phyloseq::tax_table(cbind(
#'     data_fungi_mini@tax_table,
#'     "taxon" = taxa_names(data_fungi_mini)
#'   ))
#'
#'   res_time <- ancombc_pq(
#'     data_fungi_mini,
#'     fact = "Time",
#'     levels_fact = c("0", "15"),
#'     tax_level = "taxon",
#'     verbose = TRUE
#'   )
#'
#'   signif_ancombc(res_time)
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `ANCOMBC::ancombc2()` if you
#'   use this function.
signif_ancombc <- function(ancombc_res,
                           filter_passed = TRUE,
                           filter_diff = TRUE,
                           min_abs_lfc = 0) {
  signif_ancombc_res <- ancombc_res$res
  clnames <- colnames(signif_ancombc_res)

  name_modality <-
    gsub(
      "passed_ss", "",
      clnames[grepl("passed_ss", clnames) &
        !grepl("Intercept", clnames)]
    )

  if (filter_passed) {
    signif_ancombc_res <- signif_ancombc_res %>%
      filter(.data[[paste0("passed_ss", name_modality)]])
  }

  if (filter_diff) {
    signif_ancombc_res <- signif_ancombc_res %>%
      filter(.data[[paste0("diff", name_modality)]])
  }
  signif_ancombc_res <- signif_ancombc_res %>%
    filter(abs(.data[[paste0("lfc", name_modality)]]) > min_abs_lfc)

  return(signif_ancombc_res)
}
################################################################################


################################################################################
#' Plot ANCOMBC2 result for phyloseq object
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams clean_pq
#' @param ancombc_res (required) the result of the ancombc_pq function
#'   For the moment only bimodal factors are possible.
#' @param filter_passed (logical, default TRUE) Do we filter using the column
#'   passed_ss? The passed_ss value is  TRUE if the taxon passed the sensitivity
#'   analysis, i.e., adding different pseudo-counts to 0s would not change the results.
#' @param filter_diff (logical, default TRUE) Do we filter using the column
#'   diff? The diff value is TRUE if the taxon is significant
#'   (has q less than alpha)
#' @param min_abs_lfc (integer, default 0) Minimum absolute value to filter
#'   results based on Log Fold Change. For ex. a value of 1 filter out taxa
#'   for which the abundance in a given level of the modality is not at least
#'   the double of the abundance in the other level.
#' @param tax_col The taxonomic level (must be present in `tax_table` slot)
#'   to color the points
#' @param tax_label The taxonomic level (must be present in `tax_table` slot)
#'   to add label
#' @param add_marginal_vioplot (logical, default TRUE) Do we add a marginal
#'   vioplot representing all the taxa lfc from ancombc_res.
#' @param add_label (logical, default TRUE) Do we add a label?
#' @param add_hline_cut_lfc (logical, default NULL) Do we add two horizontal
#'   lines when min_abs_lfc is set (different from zero)?
#'
#' @return A ggplot2 object. If add_marginal_vioplot is TRUE, this is a
#'   patchworks of plot made using `patchwork::plot_layout()`.
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("mia")) {
#'   data_fungi_mini@tax_table <- phyloseq::tax_table(cbind(
#'     data_fungi_mini@tax_table,
#'     "taxon" = taxa_names(data_fungi_mini)
#'   ))
#'
#'   res_time <- ancombc_pq(
#'     data_fungi_mini,
#'     fact = "Time",
#'     levels_fact = c("0", "15"),
#'     tax_level = "taxon",
#'     verbose = TRUE
#'   )
#'
#'   plot_ancombc_pq(data_fungi_mini, res_time,
#'     filter_passed = FALSE,
#'     tax_label = "Genus", tax_col = "Order"
#'   )
#'   plot_ancombc_pq(data_fungi_mini, res_time, tax_col = "Genus")
#'   plot_ancombc_pq(data_fungi_mini, res_time,
#'     filter_passed = FALSE,
#'     filter_diff = FALSE, tax_col = "Family", add_label = FALSE
#'   )
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `ANCOMBC::ancombc2()` if you
#'   use this function.
#' @author Adrien Taudière


plot_ancombc_pq <-
  function(physeq,
           ancombc_res,
           filter_passed = TRUE,
           filter_diff = TRUE,
           min_abs_lfc = 0,
           tax_col = "Genus",
           tax_label = "Species",
           add_marginal_vioplot = TRUE,
           add_label = TRUE,
           add_hline_cut_lfc = NULL) {
    verify_pq(physeq)

    if (is.null(add_hline_cut_lfc)) {
      if (min_abs_lfc != 0) {
        add_hline_cut_lfc <- TRUE
      } else {
        add_hline_cut_lfc <- FALSE
      }
    }

    signif_ancombc_res <- signif_ancombc(
      ancombc_res,
      filter_passed = filter_passed,
      filter_diff = filter_diff,
      min_abs_lfc = min_abs_lfc
    )

    clnames <- colnames(ancombc_res$res)
    name_modality <-
      gsub(
        "passed_ss", "",
        clnames[grepl("passed_ss", clnames) &
          !grepl("Intercept", clnames)]
      )

    taxtable <- data.frame(physeq@tax_table)
    taxtable$taxon <- taxa_names(physeq)

    df <-
      left_join(signif_ancombc_res, taxtable, by = join_by("taxon" == "taxon"))
    df[[tax_label]] <- gsub("unidentified", NA, df[[tax_label]])

    p <- ggplot(
      df,
      aes(
        y = reorder(taxon, .data[[paste0("lfc", name_modality)]]),
        x = .data[[paste0("lfc", name_modality)]],
        color = .data[[tax_col]]
      )
    ) +
      geom_vline(xintercept = 0) +
      geom_segment(aes(
        xend = 0,
        y = reorder(taxon, .data[[paste0("lfc", name_modality)]]),
        yend = reorder(taxon, .data[[paste0("lfc", name_modality)]])
      ), color = "darkgrey") +
      geom_point(size = 3)

    if (add_label) {
      p <- p + geom_label(aes(
        label = .data[[tax_label]], x =
          0
      ))
    }
    if (add_hline_cut_lfc) {
      p <- p +
        geom_vline(
          xintercept = min_abs_lfc,
          color = "grey20",
          lty = 2
        ) +
        geom_vline(
          xintercept = -min_abs_lfc,
          color = "grey20",
          lty = 2
        )
    }

    if (add_marginal_vioplot) {
      marg_vio <-
        (
          ggplot(data = ancombc_res$res, aes(x = .data[[paste0("lfc", name_modality)]], y = 1)) +
            geom_violin(fill = transp("grey", 0.5)) +
            geom_point(aes(color = .data[[paste0("diff", name_modality)]] &
              .data[[paste0("passed_ss", name_modality)]]), position = position_jitter()) +
            scale_color_manual(values = c("grey70", "grey20")) +
            theme(legend.position = "none")
        )

      p <-
        marg_vio / (p + xlim(layer_scales(marg_vio)$x$get_limits())) +
        patchwork::plot_layout(heights = c(1, 7), axes = "collect")
    }
    return(p)
  }
################################################################################

################################################################################
#' Show taxa which are present in only one given level of a modality
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams clean_pq
#' @param modality (required) The name of a column present in the `@sam_data` slot
#'   of the physeq object. Must be a character vector or a factor.
#' @param level (required) The level (must be present in modality) of interest
#' @param min_nb_seq_taxa (default 0 = no filter) The minimum number of sequences per taxa
#' @param min_nb_samples_taxa (default 0 = no filter) The minimum number of samples per taxa
#'
#' @return A vector of taxa names
#' @export
#'
#' @author Adrien Taudière
#' @examples
#' data_fungi_mini_woNA4height <- subset_samples(
#'   data_fungi_mini,
#'   !is.na(data_fungi_mini@sam_data$Height)
#' )
#' taxa_only_in_one_level(data_fungi_mini_woNA4height, "Height", "High")
taxa_only_in_one_level <- function(physeq,
                                   modality,
                                   level,
                                   min_nb_seq_taxa = 0,
                                   min_nb_samples_taxa = 0) {
  if (min_nb_seq_taxa > 0) {
    physeq <-
      subset_taxa_pq(physeq, taxa_sums(physeq) >= min_nb_seq_taxa)
  }
  if (min_nb_samples_taxa > 0) {
    physeq <-
      subset_taxa_pq(
        physeq,
        taxa_sums(as_binary_otu_table(physeq)) >= min_nb_samples_taxa
      )
  }

  physeq_merged <- clean_pq(merge_samples2(physeq, modality))

  physeq_merged_only_one_level <-
    subset_taxa_pq(physeq_merged, taxa_sums(as_binary_otu_table(physeq_merged)) ==
      1)
  physeq_merged_only_level_given <-
    clean_pq(subset_samples_pq(
      physeq_merged_only_one_level,
      rownames(physeq_merged_only_one_level@sam_data) == level
    ))
  return(taxa_names(physeq_merged_only_level_given))
}
################################################################################

################################################################################
#' Distribution of sequences across a factor for one taxon
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot
#'    different lines
#' @param taxa_name (required): the name of the taxa
#' @param digits (default = 2) integer indicating the number of decimal places
#'   to be used (see `?round` for more information)
#'
#' @return a dataframe with levels as rows and information as column :
#'   - the number of sequences of the taxa (nb_seq)
#'   - the number of samples of the taxa (nb_samp)
#'   - the mean (mean_nb_seq) and standard deviation (sd_nb_seq) of the *nb_seq*
#'   - the mean (mean_nb_seq_when_present) *nb_seq* excluding samples with zero
#'   - the total number of samples (nb_total_samp)
#'   - the proportion of samples with the taxa
#'
#' @export
#' @author Adrien Taudière
#' @examples
#' distri_1_taxa(data_fungi, "Height", "ASV2")
#' distri_1_taxa(data_fungi, "Time", "ASV81", digits = 1)
#' @importFrom stats sd
distri_1_taxa <- function(physeq, fact, taxa_name, digits = 2) {
  physeq <- taxa_as_rows(physeq)
  df <-
    data.frame(
      "nb_seq" = tapply(
        as.vector(physeq@otu_table[taxa_name, ]),
        physeq@sam_data[[fact]], sum
      ),
      "nb_samp" = tapply(
        as.vector(as_binary_otu_table(physeq)@otu_table[taxa_name, ]),
        physeq@sam_data[[fact]],
        sum
      ),
      "mean_nb_seq" = round(tapply(
        as.vector(physeq@otu_table[taxa_name, ]),
        physeq@sam_data[[fact]], mean
      ), digits = digits),
      "sd_nb_seq" = round(tapply(
        as.vector(physeq@otu_table[taxa_name, ]),
        physeq@sam_data[[fact]], sd
      ), digits = digits)
    ) %>%
    mutate("mean_nb_seq_when_present" = round(nb_seq / nb_samp,
      digits =
        digits
    ))

  df$nb_total_samp <- table(physeq@sam_data[[fact]])
  df$prop_samp <-
    round(df$nb_samp / df$nb_total_samp, digits = digits)
  return(df)
}
################################################################################

################################################################################
#' Partition the Variation of a phyloseq object by 2, 3, or 4 Explanatory Matrices
#' @description
#' `r lifecycle::badge("experimental")`
#'   The function partitions the variation in otu_table using
#'   distance (Bray per default) with respect to two, three, or four explanatory
#'   tables, using
#'   adjusted R² in redundancy analysis ordination (RDA) or distance-based
#'   redundancy analysis. If response is a single vector, partitioning is by
#'   partial regression. Collinear variables in the explanatory tables do NOT
#'   have to be removed prior to partitioning. See [vegan::varpart()] for more
#'   information.
#'
#' @inheritParams clean_pq
#' @param list_component (required) A named list of 2, 3 or four vectors with
#'   names from the `@sam_data` slot.
#' @param dist_method (default "bray") the distance used. See
#'   [phyloseq::distance()] for all available distances or run
#'   [phyloseq::distanceMethodList()].
#'   For "aitchison" and "robust.aitchison" distance, [vegan::vegdist()]
#'   function is directly used.
#' @param dbrda_computation (logical) Do dbrda computations are runned for each
#'  individual component (each name of the list component) ?
#'
#' @return an object of class "varpart", see [vegan::varpart()]
#' @export
#' @author Adrien Taudière
#' @seealso [var_par_rarperm_pq()], [vegan::varpart()], [plot_var_part_pq()]
#' @examples
#' \donttest{
#' if (requireNamespace("vegan")) {
#'   data_fungi_woNA <-
#'     subset_samples(data_fungi, !is.na(Time) & !is.na(Height))
#'   res_var <- var_par_pq(data_fungi_woNA,
#'     list_component = list(
#'       "Time" = c("Time"),
#'       "Size" = c("Height", "Diameter")
#'     ),
#'     dbrda_computation = TRUE
#'   )
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `vegan::varpart()` if you
#'   use this function.
var_par_pq <-
  function(physeq,
           list_component,
           dist_method = "bray",
           dbrda_computation = TRUE) {
    physeq <- taxa_as_columns(physeq)

    verify_pq(physeq)

    if (dist_method %in% c("robust.aitchison", "aitchison")) {
      dist_physeq <-
        vegan::vegdist(as(otu_table(physeq, taxa_are_rows = FALSE), "matrix"),
          method = dist_method
        )
    } else {
      dist_physeq <- phyloseq::distance(physeq, method = dist_method)
    }

    for (i in 1:length(list_component)) {
      assign(
        names(list_component)[i],
        as.data.frame(unclass(physeq@sam_data[, list_component[[i]]]))
      )
    }

    if (length(list_component) == 2) {
      res_varpart <-
        vegan::varpart(
          dist_physeq,
          eval(sym(names(list_component)[1])),
          eval(sym(names(list_component)[2]))
        )
    } else if (length(list_component) == 3) {
      res_varpart <-
        vegan::varpart(
          dist_physeq,
          eval(sym(names(list_component)[1])),
          eval(sym(names(list_component)[2])),
          eval(sym(names(list_component)[3]))
        )
    } else if (length(list_component) == 4) {
      res_varpart <-
        vegan::varpart(
          dist_physeq,
          eval(sym(names(list_component)[1])),
          eval(sym(names(list_component)[2])),
          eval(sym(names(list_component)[3])),
          eval(sym(names(list_component)[4]))
        )
    } else {
      stop("The list_component must be of length 2, 3 or 4")
    }

    if (dbrda_computation) {
      res_varpart$dbrda_result <- list()
      for (i in 1:length(list_component)) {
        res_varpart$dbrda_result[[i]] <-
          anova(vegan::dbrda(
            as.formula(paste0(
              "dist_physeq ~ ", paste(c(list_component[[i]]), collapse = " + ")
            )),
            data = as.data.frame(unclass(
              physeq@sam_data
            ))
          ))
      }
    }
    res_varpart$Xnames <- names(list_component)
    return(res_varpart)
  }
################################################################################


################################################################################
#' Partition the Variation of a phyloseq object with rarefaction permutations
#' @description
#' `r lifecycle::badge("experimental")`
#'
#'   This is an extension of the function [var_par_pq()]. The main addition is
#'   the computation of nperm permutations with rarefaction even depth by
#'   sample. The return object
#'
#'
#' @inheritParams clean_pq
#' @param list_component (required) A named list of 2, 3 or four vectors with
#'   names from the `@sam_data` slot.
#' @param dist_method (default "bray") the distance used. See
#'   [phyloseq::distance()] for all available distances or run
#'   [phyloseq::distanceMethodList()].
#'   For aitchison and robust.aitchison distance, [vegan::vegdist()]
#'   function is directly used.#' @param fill_bg
#' @param nperm (int) The number of permutations to perform.
#' @param quantile_prob (float, `[0:1]`) the value to compute the quantile.
#'   Minimum quantile is compute using 1-quantile_prob.
#' @param dbrda_computation (logical) Do dbrda computations are runned for each
#'  individual component (each name of the list component) ?
#' @param dbrda_signif_pval (float, `[0:1]`) The value under which the dbrda is
#'   considered significant.
#' @param sample.size (int) A single integer value equal to the number of
#'   reads being simulated, also known as the depth. See
#'   [phyloseq::rarefy_even_depth()].
#' @param verbose (logical). If TRUE, print additional informations.
#' @param progress_bar (logical, default TRUE) Do we print progress during
#'   the calculation?
#'
#' @return A list of class varpart with additional information in the
#'  `$part$indfract` part. Adj.R.square is the mean across permutation.
#'   Adj.R.squared_quantil_min and Adj.R.squared_quantil_max represent
#'   the quantile values of adjuste R squared
#' @export
#' @seealso [var_par_pq()], [vegan::varpart()], [plot_var_part_pq()]
#' @author Adrien Taudière
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
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `vegan::varpart()` if you
#'   use this function.
var_par_rarperm_pq <-
  function(physeq,
           list_component,
           dist_method = "bray",
           nperm = 99,
           quantile_prob = 0.975,
           dbrda_computation = FALSE,
           dbrda_signif_pval = 0.05,
           sample.size = min(sample_sums(physeq)),
           verbose = FALSE,
           progress_bar = TRUE) {
    physeq <- taxa_as_columns(physeq)
    verify_pq(physeq)

    if (progress_bar) {
      pb <- txtProgressBar(
        min = 0,
        max = nperm,
        style = 3,
        width = 50,
        char = "="
      )
    }

    if (dist_method %in% c("robust.aitchison", "aitchison")) {
      dist_physeq <-
        vegdist(as(otu_table(physeq, taxa_are_rows = FALSE), "matrix"), method = dist_method)
    } else {
      dist_physeq <- phyloseq::distance(physeq, method = dist_method)
    }

    res_perm <- list()
    for (i in 1:nperm) {
      res_perm[[i]] <-
        var_par_pq(
          physeq =
            rarefy_even_depth(
              physeq,
              rngseed = i,
              sample.size = sample.size,
              verbose = verbose
            ),
          list_component = list_component,
          dist_method = dist_method,
          dbrda_computation = dbrda_computation
        )

      if (progress_bar) {
        setTxtProgressBar(pb, i)
      }
    }
    res_varpart <- var_par_pq(
      physeq = physeq,
      list_component = list_component,
      dist_method = dist_method,
      dbrda_computation = dbrda_computation
    )

    if (dbrda_computation) {
      res_varpart$dbrda_result_prop_pval_signif <-
        rowSums(sapply(res_perm, function(x) {
          sapply(x$dbrda_result, function(xx) {
            xx$`Pr(>F)`[[1]]
          })
        }) < dbrda_signif_pval) / nperm
    }


    res_varpart$part$indfract$R.square <-
      rowMeans(sapply(res_perm, function(x) {
        (x$part$indfract$R.square)
      }))
    res_varpart$part$indfract$R.square_quantil_max <-
      apply(sapply(res_perm, function(x) {
        (x$part$indfract$R.square)
      }), 1, function(xx) {
        quantile(xx, probs = quantile_prob, na.rm = TRUE)
      })
    res_varpart$part$indfract$R.square_quantil_min <-
      apply(sapply(res_perm, function(x) {
        (x$part$indfract$R.square)
      }), 1, function(xx) {
        quantile(xx, probs = 1 - quantile_prob, na.rm = TRUE)
      })

    res_varpart$part$indfract$Adj.R.square <-
      rowMeans(sapply(res_perm, function(x) {
        (x$part$indfract$Adj.R.square)
      }))
    res_varpart$part$indfract$Adj.R.squared_quantil_max <-
      apply(sapply(res_perm, function(x) {
        (x$part$indfract$Adj.R.square)
      }), 1, function(xx) {
        quantile(xx, probs = quantile_prob, na.rm = TRUE)
      })
    res_varpart$part$indfract$Adj.R.squared_quantil_min <-
      apply(sapply(res_perm, function(x) {
        (x$part$indfract$Adj.R.square)
      }), 1, function(xx) {
        quantile(xx, probs = 1 - quantile_prob, na.rm = TRUE)
      })
    return(res_varpart)


    res_varpart$Xnames <- names(list_component)
    return(res_varpart)
  }
################################################################################
