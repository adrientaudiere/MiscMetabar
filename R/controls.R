



################################################################################
#' Calculate ecological distance among positive controls vs
#'   distance for all samples
#' @aliases dist_pos_control
#' @details  Compute distance among positive controls,
#'   i.e. samples which are duplicated
#'   to test for variation, for example in
#'    (i) a step in the sampling,
#'    (ii) a step in the extraction,
#'    (iii) a step in the sequencing.
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param samples_names (required): a vector of names for samples with
#'   positives controls of the same samples having the same name
#' @param method (default = "bray"): a method to calculate
#'    the distance, parsed to vegdist
#' @return A list of two dataframes with
#'   (i) the distance among positive controls and
#'   (ii) the distance among all samples
#' @export
#'
#' @author Adrien Taudière

dist_pos_control <- function(physeq, samples_names, method = "bray") {

  res <- list()
  dist_control <- vector()

  for (i in levels(samples_names)) {
    interm <- physeq@otu_table[samples_names ==  i, ]
    if (dim(interm)[1] > 1) {
      dist_control <-
        c(dist_control, as.numeric(vegan::vegdist(interm, method = method)))
    }
    else {
      dist_control <- c(dist_control, NA)
    }
  }
  names(dist_control) <-  levels(samples_names)
  res$dist_control_samples <- data.frame("dist_control" =
                                         na.omit(dist_control))

  # Compute distance among all samples
  if (taxa_are_rows(physeq)) {
    matdist_interm <- t(physeq@otu_table)
  }
  else {
    matdist_interm <- physeq@otu_table
  }

  res[[1]] <-
    data.frame("distAllSamples" = as.vector(vegdist(
      matdist_interm, diag = FALSE, upper = TRUE
    )))
  res[[2]] <-
    data.frame("dist_control_samples" = as.vector(na.omit(dist_control)))
  names(res) <- c("distAllSamples", "dist_controlontrolSamp
les")

  return(res)
}
################################################################################



################################################################################
#' Subset taxa using a taxa control (e.g. truffle root tips) through 3 methods.
#' @aliases subset_taxa_tax_control
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param taxa_distri (required): a vector of length equal to the number of
#'   samples with the number of sequences per samples for the taxa control
#' @param method (default = "mean"): a method to calculate the cut-off value.
#'   There is 6 methods:
#'   (i)   cutoff_seq: discard taxa with less than the number of sequence
#'           than taxa control,
#'   (ii)  cutoff_mixt: using mixture models,
#'   (iii) cutoff_diff: using a minimum difference threshold
#'           (need the argument min_diff_for_cutoff)
#'   (iv)  min: the minimum of the three firsts methods
#'   (v)   max: the maximum of the three firsts methods
#'   (vi)  mean: the mean of the three firsts methods
#' @param min_diff_for_cutoff (default = 10): argument for method cutoff_diff
#' @return A new \code{\link{phyloseq-class}} object.
#' @export
#'
#' @author Adrien Taudière
subset_taxa_tax_control <-
  function(physeq,
           taxa_distri,
           method = "mean",
           min_diff_for_cutoff = 10) {
    cutoff_seq <- vector(mode = "numeric", length = nsamples(physeq))
    cutoff_mixt <-
      vector(mode = "numeric", length = nsamples(physeq))
    cutoff_diff <-
      vector(mode = "numeric", length = nsamples(physeq))
    cutoffs <- vector(mode = "numeric", length = nsamples(physeq))


    for (i in 1:nsamples(physeq)) {
      # for each samples

      if (method %in% c("min", "max", "mean", "cutoff_mixt")) {
        find_cutoff <- function(proba = 0.5, il = index_lower) {
          ## Cutoff such that Pr[drawn from bad component] == proba
          f <- function(x) {
            proba - (
              model$lambda[il] * dnorm(x, model$mu[il], model$sigma[il]) /
                (
                  model$lambda[1] * dnorm(x, model$mu[1], model$sigma[1]) +
                    model$lambda[2] * dnorm(x, model$mu[2], model$sigma[2])
                )
            )
          }
          return(uniroot(
            f = f,
            lower = 1,
            upper = 1000
          )$root)  # Careful with division by zero if changing lower and upper
        }

        physeq_mixture <- as.vector(physeq@otu_table[i, ])
        x_mixture <- physeq_mixture[physeq_mixture > 0]
        if (length(x_mixture) > 25) {
          try(model <-
                mixtools::normalmixEM(
                  x = x_mixture,
                  k = 2,
                  epsilon = 1e-03,
                  maxrestarts = 1000
                ),
              silent = TRUE)
          try(index_lower <-
                which.min(model$mu), silent = TRUE)
          # Index of component with lower mean)

          cutoff_mixt[i] <-
            try(find_cutoff(proba = 0.5,  il = index_lower))
        } else {
          cutoff_mixt[i] <- NA
        }
        cutoff_mixt <- as.numeric(cutoff_mixt)
      }

      if (method %in% c("min", "max", "mean", "cutoff_diff")) {
        physeq_diff <- sort(as.vector(as.vector(physeq@otu_table[i, ])))
        cond <- diff(physeq_diff) > min_diff_for_cutoff
        cutoff_diff[i] <- min(physeq_diff[cond])
      }

      if (method %in% c("min", "max", "mean", "cutoff_seq")) {
        physeq_seq <- sort(as.vector(as.vector(physeq@otu_table[i, ])))
        cutoff_seq[i] <- min(physeq_seq[physeq_seq > taxa_distri[i]])
        if (is.infinite(cutoff_seq[i])) {
          cutoff_seq[i] <- taxa_distri[i]
        }
      }
    }

    if (method == "cutoff_mixt")  {
      cutoffs <- cutoff_mixt
    } else if (method == "cutoff_diff")  {
      cutoffs <- cutoff_diff
    } else  if (method == "cutoff_seq")  {
      cutoffs <-  cutoff_seq
    } else  if (method == "mean")  {
      cutoffs <-
        apply(cbind(cutoff_mixt, cutoff_diff, cutoff_seq), 1, mean, na.rm = T)
    } else  if (method == "min")  {
      cutoffs <-
        apply(cbind(cutoff_mixt, cutoff_diff, cutoff_seq), 1, min, na.rm = T)
    } else  if (method == "max")  {
      cutoffs <-
        apply(cbind(cutoff_mixt, cutoff_diff, cutoff_seq), 1, max, na.rm = T)
    } else {
      stop("The method name is not valid.")
    }

    new_physeq <- physeq
    for (i in 1:nsamples(new_physeq)) {
      new_physeq@otu_table[i, new_physeq@otu_table[i, ] < floor(cutoffs)[i]] <-
        0
    }

    new_physeq <-
      prune_taxa(colSums(new_physeq@otu_table) > 0, new_physeq)
    ntaxa_discard <- ntaxa(physeq) - ntaxa(new_physeq)
    nseq_discard <-
      sum(physeq@otu_table) - sum(new_physeq@otu_table)
    message(
      paste(
        "The filtering processes discard",
        ntaxa_discard,
        "taxa and",
        nseq_discard,
        "sequences",
        sep = " "
      )
    )
    return(new_physeq)
  }
