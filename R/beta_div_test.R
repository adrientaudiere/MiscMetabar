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
#'   which samples to merge using the
#'   \code{\link[speedyseq]{merge_samples2}} function.
#'   Need to be in \code{physeq@sam_data}
#' @param nperm (int) The number of permutations to perform.
#' @param return_plot (logical) Do we return only the result
#'   of the test or do we plot the result?
#' @param title The title of the Graph.
#' @param ... other params for be passed on
#'   [phyloseqGraphTest::graph_perm_test()] function
#'
#' @examples
#' data(enterotype)
#' graph_test_pq(enterotype, fact = "SeqTech")
#'
#' clean_enterotype <- subset_samples(
#'   enterotype,
#'   !is.na(enterotype@sam_data$Enterotype)
#' )
#' graph_test_pq(clean_enterotype, fact = "Enterotype")
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
    physeq <- clean_pq(physeq)
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


#' @title Permanova on a phyloseq object
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A wrapper for the [vegan::adonis2()] function in the case of `physeq` object.
#' @inheritParams clean_pq
#' @param formula (required) the right part of a formula for [vegan::adonis2()].
#'   Variables must be present in the `physeq@sam_data` slot.
#' @param dist_method (default "bray") the distance used. See 
#'   [phyloseq:::distance()] for all available distances or run
#'   `phyloseq::distanceMethodList()`.
#'   For aitchison and robust.aitchison distance, [vegan::vegdist()] 
#'   function is directly used.
#' @param merge_sample_by a vector to determine
#'   which samples to merge using the [speedyseq::merge_samples2()]
#'   function. Need to be in `physeq@sam_data`
#' @param na_remove (logical, default FALSE) If set to TRUE, remove samples with 
#'   NA in the variables set in formula.
#' @param correction_for_sample_size (logical, default FALSE) If set to TRUE, the 
#'   sample size (number of sequences by samples) is add to formula in the form 
#'   `y~Library_Size + Biological_Effect` following recommendation of 
#'   [Weiss et al. 2017](https://doi.org/10.1186/s40168-017-0237-y). 
#'   `correction_for_sample_size` overcome `rarefy_nb_seqs` if both are TRUE.
#' @param rarefy_nb_seqs (logical, default FALSE) Rarefy each sample
#'   (before merging if merge_sample_by is set) using `phyloseq::rarefy_even_depth()`.
#'   if `correction_for_sample_size` is TRUE, rarefy_nb_seqs will have no effect.
#' @return The function returns an anova.cca result object with a
#'   new column for partial R^2. See help of [vegan::adonis2()] for more information.
#'
#' @examples
#' data(enterotype)
#' adonis_pq(enterotype, "SeqTech*Enterotype", na_remove = TRUE)
#' adonis_pq(enterotype, "SeqTech")
#' adonis_pq(enterotype, "SeqTech", dist_method = "jaccard")
#' adonis_pq(enterotype, "SeqTech", dist_method = "robust.aitchison")
#' @export
#' @importFrom stats reformulate
#' @author Adrien Taudière

adonis_pq <- function(physeq,
  formula,
  dist_method = "bray",
  merge_sample_by = NULL,
  na_remove = FALSE,
  correction_for_sample_size = FALSE,
  rarefy_nb_seqs = FALSE) {

  physeq <- clean_pq(
      physeq,
      force_taxa_as_columns = TRUE,
      remove_empty_samples = TRUE,
      remove_empty_taxa = FALSE,
      clean_samples_names = FALSE,
      silent = TRUE
    )

  if(dist_method %in% c("aitchison", "robust.aitchison")){
    phy_dist <- paste0('vegan::vegdist(as.matrix(physeq@otu_table), method="', dist_method, '")')
  } else {
    phy_dist <- paste0('phyloseq:::distance(physeq, method="', dist_method, '")')
  }
  
  .formula <- reformulate(formula, response = phy_dist)
  termf <- terms(.formula)
  term_lab <- attr(termf, "term.labels")[attr(termf, "order") == 1]

  verify_pq(physeq)

  if(na_remove){
    new_physeq <- physeq
    for(tl in term_lab) {
      new_physeq <- subset_samples(new_physeq, !is.na(physeq@sam_data[[tl]]))
    }
    if(nsamples(physeq)-nsamples(new_physeq) > 0){
      message(paste0(nsamples(physeq)-nsamples(new_physeq), 
                   " were discarded due to NA in variables present in formula."))
    }
    physeq <- new_physeq
  }

  if (!is.null(merge_sample_by)) {
    physeq <- speedyseq::merge_samples2(physeq, merge_sample_by)
    physeq <- clean_pq(physeq)
  }    

  if (correction_for_sample_size) {
    formula <- paste0("sample_size+", formula)
    .formula <- reformulate(formula, response = phy_dist)
  } else if (rarefy_nb_seqs) {
    physeq <- rarefy_even_depth(physeq)
    physeq <- clean_pq(physeq)
  }

  verify_pq(physeq)
  metadata <- as(sample_data(physeq), "data.frame")
  if (correction_for_sample_size) { 
    metadata$sample_size <- sample_sums(physeq)
  }

  res_ado <- vegan::adonis2(.formula, data = metadata)
  return(res_ado)
}
