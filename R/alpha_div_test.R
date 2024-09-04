################################################################################
#' Calculate hill number and compute Tuckey post-hoc test
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#'
#' Note that, by default, this function use a sqrt of the read numbers in the linear
#'   model in order to correct for uneven sampling depth.
#' @aliases hill_tuckey_pq
#' @inheritParams clean_pq
#' @param modality (required) the variable to test
#' @param hill_scales (a vector of integer) The list of q values to compute
#'   the hill number H^q. If Null, no hill number are computed. Default value
#'   compute the Hill number 0 (Species richness), the Hill number 1
#'   (exponential of Shannon Index) and the Hill number 2 (inverse of Simpson
#'   Index).
#' @param silent (logical) If TRUE, no message are printing.
#' @param correction_for_sample_size (logical, default TRUE) This function
#'   use a sqrt of the read numbers in the linear model in order to
#'   correct for uneven sampling depth.
#' @return A ggplot2 object
#'
#' @export
#'
#' @author Adrien Taudière
#' @examples
#' data("GlobalPatterns", package = "phyloseq")
#' GlobalPatterns@sam_data[, "Soil_logical"] <-
#'   ifelse(GlobalPatterns@sam_data[, "SampleType"] == "Soil", "Soil", "Not Soil")
#' hill_tuckey_pq(GlobalPatterns, "Soil_logical")
#' hill_tuckey_pq(GlobalPatterns, "Soil_logical", hill_scales = 1:2)
hill_tuckey_pq <- function(
    physeq,
    modality,
    hill_scales = c(0, 1, 2),
    silent = TRUE,
    correction_for_sample_size = TRUE) {
  modality_vector <-
    as.factor(as.vector(unlist(unclass(physeq@sam_data[, modality]))))

  if (length(modality_vector) != dim(physeq@otu_table)[2]) {
    physeq@otu_table <- t(physeq@otu_table)
  }
  read_numbers <- apply(physeq@otu_table, 2, sum)

  physeq <- taxa_as_rows(physeq)
  otu_hill <-
    vegan::renyi(t(physeq@otu_table),
      scales = hill_scales,
      hill = TRUE
    )

  colnames(otu_hill) <- paste0("Hill_", hill_scales)
  tuk <- list()
  for (i in seq_along(hill_scales)) {
    if (correction_for_sample_size) {
      tuk[[i]] <-
        stats::TukeyHSD(stats::aov(lm(otu_hill[, i] ~ sqrt(read_numbers))$residuals ~ modality_vector))
    } else {
      tuk[[i]] <-
        stats::TukeyHSD(stats::aov(otu_hill[, i] ~ modality_vector))
    }
  }
  df <- do.call(
    "rbind",
    sapply(tuk, function(x) {
      data.frame(x$modality_vector)
    }, simplify = FALSE)
  )
  colnames(df) <- colnames(tuk[[1]]$modality_vector)
  df$x <- paste0(
    "Hill_",
    c(
      sort(rep(hill_scales, dim(
        tuk[[1]]$modality_vector
      )[1]))
    ), "__",
    rownames(tuk[[1]]$modality_vector)
  )

  df$modality <- rownames(tuk[[1]]$modality_vector)

  p <- ggplot(data = df) +
    geom_linerange(aes(ymax = upr, ymin = lwr, x = x), linewidth = 2) +
    geom_point(aes(x = x, y = diff),
      size = 4,
      shape = 21,
      fill = "white"
    ) +
    coord_flip() +
    theme_gray() +
    geom_hline(yintercept = 0) +
    ylab("Differences in mean levels (value and confidence intervals at 95%)") +
    xlab("") +
    ggtitle("Results of the Tuckey HSD testing for differences
    in mean Hill numbers")

  return(p)
}
################################################################################


################################################################################
#' Test multiple times effect of factor on Hill diversity
#'   with different rarefaction even depth
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' This reduce the risk of a random drawing of a exceptional situation of an unique rarefaction.
#' @inheritParams clean_pq
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot
#'    different lines
#' @param hill_scales (a vector of integer) The list of q values to compute
#'   the hill number H^q. If Null, no hill number are computed. Default value
#'   compute the Hill number 0 (Species richness), the Hill number 1
#'   (exponential of Shannon Index) and the Hill number 2 (inverse of Simpson
#'   Index).
#' @param nperm (int) The number of permutations to perform.
#' @param sample.size (int) A single integer value equal to the number of
#'   reads being simulated, also known as the depth. See
#'   [phyloseq::rarefy_even_depth()].
#' @param verbose (logical). If TRUE, print additional informations.
#' @param progress_bar (logical, default TRUE) Do we print progress during
#'   the calculation?
#' @param p_val_signif (float, `[0:1]`) The mimimum value of p-value to count a
#'   test as significant int the `prop_signif` result.
#' @param type A character specifying the type of statistical approach
#'   (See [ggstatsplot::ggbetweenstats()] for more details):
#'
#'   - "parametric"
#'   - "nonparametric"
#'   - "robust"
#'   - "bayes"
#'
#' @param ... Other arguments passed on to [ggstatsplot::ggbetweenstats()] function
#' @seealso [ggstatsplot::ggbetweenstats()], [hill_pq()]
#' @return A list of 6 components :
#'
#' - method
#' - expressions
#' - plots
#' - pvals
#' - prop_signif
#' - statistics
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggstatsplot")) {
#'   hill_test_rarperm_pq(data_fungi, "Time", nperm = 2)
#'   res <- hill_test_rarperm_pq(data_fungi, "Height", nperm = 9, p.val = 0.9)
#'   patchwork::wrap_plots(res$plots[[1]])
#'   res$plots[[1]][[1]] + res$plots[[2]][[1]] + res$plots[[3]][[1]]
#'   res$prop_signif
#'   res_para <- hill_test_rarperm_pq(data_fungi, "Height", nperm = 9, type = "parametrique")
#'   res_para$plots[[1]][[1]] + res_para$plots[[2]][[1]] + res_para$plots[[3]][[1]]
#'   res_para$pvals
#'   res_para$method
#'   res_para$expressions[[1]]
#' }
#' }
hill_test_rarperm_pq <- function(physeq,
                                 fact,
                                 hill_scales = c(0, 1, 2),
                                 nperm = 99,
                                 sample.size = min(sample_sums(physeq)),
                                 verbose = FALSE,
                                 progress_bar = TRUE,
                                 p_val_signif = 0.05,
                                 type = "non-parametrique",
                                 ...) {
  verify_pq(physeq)
  res_perm <- list()
  p_perm <- list()
  if (progress_bar) {
    pb <- txtProgressBar(
      min = 0,
      max = nperm * length(hill_scales),
      style = 3,
      width = 50,
      char = "="
    )
  }
  for (i in 1:nperm) {
    if (verbose) {
      psm <-
        psmelt_samples_pq(
          physeq = rarefy_even_depth(
            physeq,
            rngseed = i,
            sample.size = sample.size,
            verbose = verbose
          ),
          hill_scales = hill_scales
        )
    } else {
      psm <-
        suppressMessages(psmelt_samples_pq(
          physeq = rarefy_even_depth(
            physeq,
            rngseed = i,
            sample.size = sample.size,
            verbose = verbose
          ),
          hill_scales = hill_scales
        ))
    }
    p_perm[[i]] <- list()
    res_perm[[i]] <- list()
    for (j in seq_along(hill_scales)) {
      p_perm[[i]][[j]] <-
        ggstatsplot::ggbetweenstats(psm, !!fact, !!paste0("Hill_", hill_scales[[j]]),
          type = type,
          ...
        )
      res_perm[[i]][[j]] <-
        ggstatsplot::extract_stats(p_perm[[i]][[j]])
    }
    if (progress_bar) {
      setTxtProgressBar(pb, i * length(hill_scales))
    }
  }

  method <- res_perm[[1]][[1]]$subtitle_data[, c("method", "effectsize", "conf.method")]

  expressions <- sapply(res_perm, function(x) {
    sapply(x, function(xx) {
      xx$subtitle_data$expression
    })
  })
  rownames(expressions) <- paste0("Hill_", hill_scales)
  colnames(expressions) <- paste0("ngseed", 1:nperm)

  statistics <- sapply(res_perm, function(x) {
    sapply(x, function(xx) {
      xx$subtitle_data$statistic
    })
  })
  rownames(statistics) <- paste0("Hill_", hill_scales)
  colnames(statistics) <- paste0("ngseed", 1:nperm)

  pvals <- sapply(res_perm, function(x) {
    sapply(x, function(xx) {
      xx$subtitle_data$p.value
    })
  })
  rownames(pvals) <- paste0("Hill_", hill_scales)
  colnames(pvals) <- paste0("ngseed_", 1:nperm)

  prop_signif <- rowSums(pvals < p_val_signif) / ncol(pvals)
  names(prop_signif) <- paste0("Hill_", hill_scales)
  res <-
    list(
      "method" = method,
      "expressions" = expressions,
      "plots" = p_perm,
      "pvals" = pvals,
      "prop_signif" = prop_signif,
      "statistics" = statistics
    )
  return(res)
}
################################################################################



################################################################################
#' Automated model selection and multimodel inference with (G)LMs for phyloseq
#'
#' @description
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' See [glmulti::glmulti()] for more information.
#'
#' @inheritParams clean_pq
#' @param formula (required) a formula for [glmulti::glmulti()]
#'   Variables must be present in the `physeq@sam_data` slot or be one
#'   of hill number defined in hill_scales or the variable Abundance which
#'   refer to the number of sequences per sample.
#' @param fitfunction (default "lm")
#' @param hill_scales (a vector of integer) The list of q values to compute
#'   the hill number H^q. If Null, no hill number are computed. Default value
#'   compute the Hill number 0 (Species richness), the Hill number 1
#'   (exponential of Shannon Index) and the Hill number 2 (inverse of Simpson
#'   Index).
#' @param aic_step The value between AIC scores to cut for.
#' @param confsetsize The number of models to be looked for, i.e. the size of the returned confidence set.
#' @param plotty (logical) Whether to plot the progress of the IC profile when running.
#' @param level If 1, only main effects (terms of order 1) are used to build
#'   the candidate set. If 2, pairwise interactions are also used (higher order
#'   interactions are currently ignored)
#' @param method The method to be used to explore the candidate set of models.
#'   If "h" (default) an exhaustive screening is undertaken.
#'   If "g" the genetic algorithm is employed (recommended for large candidate sets).
#'   If "l", a very fast exhaustive branch-and-bound algorithm is used.
#'   Package leaps must then be loaded, and this can only be applied to linear models
#'   with covariates and no interactions. If "d", a simple summary of the candidate set
#'   is printed, including the number of candidate models.
#' @param crit The Information Criterion to be used. Default is the small-sample corrected AIC (aicc). This should be a function that accepts a fitted model as first argument. Other provided functions are the classic AIC, the Bayes IC (bic), and QAIC/QAICc (qaic and qaicc).
#' @param ... Other arguments passed on to [glmulti::glmulti()] function
#'
#' @return A data.frame summarizing the glmulti results with columns
#'
#'  -estimates
#'  -unconditional_interval
#'  -nb_model"
#'  -importance
#'  -alpha
#' @export
#' @seealso  [glmulti::glmulti()]
#' @examples
#' \donttest{
#' if (requireNamespace("glmulti")) {
#'   res_glmulti <-
#'     glmutli_pq(data_fungi, "Hill_0 ~ Hill_1 + Abundance + Time + Height", level = 1)
#'   res_glmulti
#'   res_glmulti_interaction <-
#'     glmutli_pq(data_fungi, "Hill_0 ~ Abundance + Time + Height", level = 2)
#'   res_glmulti
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to [glmulti::glmulti()] if you
#'   use this function.
glmutli_pq <-
  function(physeq,
           formula,
           fitfunction = "lm",
           hill_scales = c(0, 1, 2),
           aic_step = 2,
           confsetsize = 100,
           plotty = FALSE,
           level = 1,
           method = "h",
           crit = "aicc",
           ...) {
    psm_samp <- psmelt_samples_pq(physeq, hill_scales = hill_scales)

    res_glmulti <- do.call(glmulti::glmulti, list(
      y = formula(formula),
      data = psm_samp,
      crit = crit,
      level = level,
      method = method,
      fitfunction = fitfunction,
      confsetsize = confsetsize,
      plotty = plotty,
      ...
    ))

    ## AICc
    top_glmulti <- glmulti::weightable(res_glmulti)
    condition_crit <- top_glmulti[[crit]] <= (min(top_glmulti[[crit]]) + aic_step)
    if (sum(condition_crit) == 0) {
      stop("None modele are selected. Try a aic_step lower or another crit")
    }
    top_glmulti <- top_glmulti[condition_crit, ]

    ## Stockage des meilleurs modèles
    cf <- data.frame(stats::coef(res_glmulti, icmethod = "Burnham"))

    colnames(cf) <-
      c(
        "estimates",
        "unconditional_interval",
        "nb_model",
        "importance",
        "alpha"
      )
    cf$variable <- rownames(cf)
    cf <- cf %>% filter(!grepl("Intercept", variable))

    if (fitfunction == "lm") {
      test <- list()
      R2__h0 <- NULL
      for (i in 1:nrow(top_glmulti)) {
        test[[i]] <- summary(res_glmulti@objects[[i]])
        R2__h0[i] <- test[[i]]$adj.r.squared
      }

      # message(paste0("Mean adjust r squared: ", round(mean(R2__h0), 3)))
    }
    return(cf)
  }
