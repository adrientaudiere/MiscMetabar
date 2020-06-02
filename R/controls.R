



################################################################################
#' Calculate ecological distance among positive controls vs
#'   distance for all samples
#' @aliases dist_pos_control
#' @notes   Compute distance among positive controls,
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
#'
#' @author Adrien Taudière

dist_pos_control <- function(physeq, samples_names, method = "bray") {

  res <- list()
  dist_control <- vector()

  for (i in levels(samples_names)) {
    interm <- physeq@otu_table[samples_names ==  i,]
    if (dim(interm)[1] > 1) {
      dist_control <-
        c(dist_control, as.numeric(vegan::vegdist(interm, method = method)))
    }
    else {
      dist_control <- c(dist_control, NA)
    }
  }
  names(dist_control) <-  levels(samples_names)
  res$dist_controlontrolSamples <- data.frame("dist_control" =
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
    data.frame("dist_controlontrolSamples" = as.vector(na.omit(dist_control)))
  names(res) <- c("distAllSamples", "dist_controlontrolSamp
les")

  return(res)
}
################################################################################



################################################################################
#' Subset taxa using a taxa control (e.g. truffle root tips) through 3 methods.
#' @aliases subset_taxa_taxControl
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param taxa_distri (required): a vector of length equal to the number of
#'   samples with the number of sequences per samples for the taxa control
#' @param method (default = "mean"): a method to calculate the cut-off value.
#'   There is 6 methods:
#'   (i)   cutoffSeq: discard taxa with less than the number of sequence
#'           than taxa control,
#'   (ii)  cutoffMixt: using mixture models,
#'   (iii) cutoffDiff: using a minimum difference threshold
#'           (need the argument min_diff_for_cutoff)
#'   (iv)  min: the minimum of the three firsts methods
#'   (v)   max: the maximum of the three firsts methods
#'   (vi)  mean: the mean of the three firsts methods
#' @param min_diff_for_cutoff (default = 10): argument for method cutoffDiff
#' @return A new \code{\link{phyloseq-class}} object.
#'
#' @author Adrien Taudière
subset_taxa_taxControl <-
  function(physeq,
           taxa_distri,
           method = "mean",
           min_diff_for_cutoff = 10) {
    cutoffSeq <- vector(mode = "numeric", length = nsamples(physeq))
    cutoffMixt <-
      vector(mode = "numeric", length = nsamples(physeq))
    cutoffDiff <-
      vector(mode = "numeric", length = nsamples(physeq))
    cutoffs <- vector(mode = "numeric", length = nsamples(physeq))


    for (i in 1:nsamples(physeq)) {
      # for each samples

      if (method %in% c("min", "max", "mean", "cutoffMixt")) {
        find.cutoff <- function(proba = 0.5, il = index.lower) {
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

        physeq_mixture <- as.vector(physeq@otu_table[i,])
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
          try(index.lower <-
                which.min(model$mu), silent = TRUE)
          # Index of component with lower mean)

          cutoffMixt[i] <-
            try(find.cutoff(proba = 0.5,  il = index.lower))
        } else {
          cutoffMixt[i] <- NA
        }
        cutoffMixt <- as.numeric(cutoffMixt)
      }

      if (method %in% c("min", "max", "mean", "cutoffDiff")) {
        physeq_Diff <- sort(as.vector(as.vector(physeq@otu_table[i,])))
        cond <- diff(physeq_Diff) > min_diff_for_cutoff
        cutoffDiff[i] <- min(physeq_Diff[cond])
      }

      if (method %in% c("min", "max", "mean", "cutoffSeq")) {
        physeq_Seq <- sort(as.vector(as.vector(physeq@otu_table[i,])))
        cutoffSeq[i] <- min(physeq_Seq[physeq_Seq > taxa_distri[i]])
        if (is.infinite(cutoffSeq[i])) {
          cutoffSeq[i] <- taxa_distri[i]
        }
      }
    }

    if (method == "cutoffMixt")  {
      cutoffs <- cutoffMixt
    } else if (method == "cutoffDiff")  {
      cutoffs <- cutoffDiff
    } else  if (method == "cutoffSeq")  {
      cutoffs <-  cutoffSeq
    } else  if (method == "mean")  {
      cutoffs <-
        apply(cbind(cutoffMixt, cutoffDiff, cutoffSeq), 1, mean, na.rm = T)
    } else  if (method == "min")  {
      cutoffs <-
        apply(cbind(cutoffMixt, cutoffDiff, cutoffSeq), 1, min, na.rm = T)
    } else  if (method == "max")  {
      cutoffs <-
        apply(cbind(cutoffMixt, cutoffDiff, cutoffSeq), 1, max, na.rm = T)
    } else {
      stop("The method name is not valid.")
    }

    new.physeq <- physeq
    for (i in 1:nsamples(new.physeq)) {
      new.physeq@otu_table[i, new.physeq@otu_table[i,] < floor(cutoffs)[i]] <-
        0
    }

    new.physeq <-
      prune_taxa(colSums(new.physeq@otu_table) > 0, new.physeq)
    ntaxa_discard <- ntaxa(physeq) - ntaxa(new.physeq)
    nseq_discard <-
      sum(physeq@otu_table) - sum(new.physeq@otu_table)
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
    return(new.physeq)
  }
