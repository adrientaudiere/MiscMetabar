
#' @title Performs graph-based permutation tests on phyloseq object
#' @description
#' #' A wrapper of \code{\link[phyloseqGraphTest]{graph_perm_test}} for quick plot with
#' important statistics
#' 
#' `r lifecycle::badge("experimental")`
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param fact (required): Name of the factor to cluster samples by modalities.
#' Need to be in \code{physeq@sam_data}. This should be a factor.
#' with two or more levels.
#' @param merge_sample_by (Default: NULL) : a vector to determine
#' which samples to merge using the \code{\link[speedyseq]{merge_samples2}} function.
#' Need to be in \code{physeq@sam_data}
#' @param  nperm (Default: 999): The number of permutations to perform.
#' @param return_plot (Default: TRUE, logical): do we return only the result
#' of the test or do we plot the result
#' @param ... other params for be passed on
#' `phyloseqGraphTest::graph_perm_test` function
#' 
#' @examples 
#' data(enterotype)
#' physeq_graph_test(enterotype, fact="SeqTech")
#' 
#' clean_enterotype <- subset_samples(enterotype, !is.na(enterotype@sam_data$Enterotype))
#' physeq_graph_test(clean_enterotype, fact="Enterotype")
#' @author Adrien Taudière
#'
#' @return a ggplot with a subtitle indicating the pvalue
#' and the number of permutations
#'
#' @export

graph_test_pq <- function(physeq,
                              fact,
                              merge_sample_by = NULL,
                              nperm = 999,
                              return_plot = TRUE,
                              title = "Graph Test",
                              ...) {
  
  verify_pq(physeq)

  if (!is.null(merge_sample_by)) {
    physeq <- speedyseq::merge_samples2(physeq, merge_sample_by)
    physeq <- clean_physeq(physeq)
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
          "pvalue = ", res_graph_test$pval,
          "(", length(res_graph_test$perm), " permutations)"
        )
      )
    return(p)
  }
}


#' Permanova on a phyloseq object
#' @description
#' `r lifecycle::badge("experimental")`
#' A wrapper for the [vegan::adonis2()] function in the case of `physeq` object. 
#' @param physeq (required): a \code{\link{phyloseq-class}} object.
#' @param formula (required): todo
#' @param merge_sample_by (Default: NULL) : a vector to determine
#' which samples to merge using the [speedyseq::merge_samples2()]
#' function. Need to be in `physeq@sam_data`
#' @return The function returns an anova.cca result object with a new column for partial R^2.
#' See ?adonis2 for more information.
#' 
#' @examples
#' data(enterotype)
#' clean_enterotype <- subset_samples(enterotype, !is.na(enterotype@sam_data$Enterotype))
#' adonis_pq(clean_enterotype, "SeqTech*Enterotype")
#' @export
#' @author Adrien Taudière

adonis_pq <- function(physeq, formula, merge_sample_by = NULL) {
  verify_pq(physeq)
  if (!is.null(merge_sample_by)) {
    physeq <- speedyseq::merge_samples2(physeq, merge_sample_by)
    physeq <- clean_physeq(physeq)
  }
  metadata <- as(sample_data(physeq), "data.frame")
  phy_dist <- 'phyloseq:::distance(physeq, method="bray")'
  .formula <- reformulate(formula, response = phy_dist)
  res_ado <- vegan::adonis2(.formula, data = metadata)
  return(res_ado)
}
