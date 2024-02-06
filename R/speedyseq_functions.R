#' Merge taxa in groups (vectorized version)
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Firstly release in the [speedyseq](https://github.com/mikemc/speedyseq/) R
#' package by Michael R. McLaren.
#'
#' Merge taxa in `x` into a smaller set of taxa defined by the vector `group`.
#' Taxa whose value in `group` is NA will be dropped. New taxa will be named
#' according to the most abundant taxon in each group (`phyloseq` and
#' `otu_table` objects) or the first taxon in each group (all other phyloseq
#' component objects).
#'
#' If `x` is a phyloseq object with a phylogenetic tree, then the new taxa will
#' be ordered as they are in the tree. Otherwise, the taxa order can be
#' controlled by the `reorder` argument, which behaves like the `reorder`
#' argument in [base::rowsum()]. `reorder = FALSE` will keep taxa in
#' the original order determined by when the member of each group first appears
#' in `taxa_names(x)`; `reorder = TRUE` will order new taxa according to their
#' corresponding value in `group`.
#'
#' The `tax_adjust` argument controls the handling of taxonomic disagreements
#' within groups. Setting `tax_adjust == 0` causes no adjustment; the taxonomy
#' of the new group is set to the archetype taxon (see below). Otherwise,
#' disagreements within a group at a given rank cause the values at lower ranks
#' to be set to `NA`. If `tax_adjust == 1` (the default), then a rank where all
#' taxa in the group are already NA is not counted as a disagreement, and lower
#' ranks may be kept if the taxa agree. This corresponds to the original
#' phyloseq behavior. If `tax_adjust == 2`, then these NAs are treated as a
#' disagreement; all ranks are set to NA after the first disagreement or NA.
#'
#' @param x A phyloseq object or component object
#' @param group A vector with one element for each taxon in `physeq` that
#' defines the new groups. see `base::rowsum()`.
#' @param reorder Logical specifying whether to reorder the taxa by their
#' `group` values. Ignored if `x` has (or is) a phylogenetic tree.
#' @param tax_adjust 0: no adjustment; 1: phyloseq-compatible adjustment; 2:
#' conservative adjustment
#' @export
#'
#' @seealso
#' Function in MiscMetabar that use this function: [asv2otu()]
#'
#' [base::rowsum()]
#'
#' [phyloseq::merge_taxa()]
#'
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X)) modified by Adrien Taudiere
setGeneric(
  "merge_taxa_vec",
  function(x,
           group,
           reorder = FALSE,
           tax_adjust = 1L) {
    standardGeneric("merge_taxa_vec")
  }
)

#' @rdname merge_taxa_vec
setMethod(
  "merge_taxa_vec", "phyloseq",
  function(x, group, reorder = FALSE, tax_adjust = 1L) {
    stopifnot(ntaxa(x) == length(group))
    stopifnot(tax_adjust %in% c(0L, 1L, 2L))
    # Warn the user if an impossible reordering is requested
    if (!is.null(x@phy_tree) & reorder) {
      warning("Can't reorder taxa if `x` has a `phy_tree`")
      reorder <- FALSE
    }
    # drop taxa with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding taxa will be dropped")
      x <- prune_taxa(!is.na(group), x)
      group <- group[!is.na(group)]
    }
    # Get the merged otu table with new taxa named by most abundant
    otu <- merge_taxa_vec(otu_table(x), group, reorder = reorder)
    # Adjust taxonomy if necessary
    if (!is.null(x@tax_table) & tax_adjust != 0) {
      tax <- merge_taxa_vec(tax_table(x), group,
        tax_adjust = tax_adjust,
        reorder = reorder
      )
      # Taxa in `tax` are in same order as in `otu` but are named by first in
      # group instead of max and so need to be renamed
      taxa_names(tax) <- taxa_names(otu)
    } else {
      tax <- NULL
    }
    # Create the new phyloseq object. Replacing the original otu_table with
    # the new, smaller table will automatically prune the taxonomy, tree, and
    # refseq to the smaller set of archetypal taxa.
    otu_table(x) <- otu
    if (!is.null(tax)) {
      tax_table(x) <- tax
    }
    return(x)
  }
)

#' @rdname merge_taxa_vec
setMethod(
  "merge_taxa_vec", "otu_table",
  function(x, group, reorder = FALSE) {
    stopifnot(ntaxa(x) == length(group))
    # Work with taxa as rows, and remember to flip back at end if needed
    needs_flip <- !taxa_are_rows(x)
    if (needs_flip) {
      x <- t(x)
    }
    # Drop taxa with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding taxa will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    # New taxa names are the most abundant taxon in each group; in the case of
    # ties, the first taxon is chosen. Original group order is maintained.
    new_names <- tibble(
      taxon = taxa_names(x),
      sum = taxa_sums(x),
      group = factor(group, levels = unique(group))
    ) %>%
      group_by(group) %>%
      mutate(archetype = taxon[which.max(sum)]) %>%
      group_by(group) %>%
      dplyr::slice_head()

    if (reorder) {
      new_names <- new_names %>% arrange(archetype)
    }
    # Compute new table with base::rowsum(). The call to rowsum() makes the
    # rownames the group names.
    otu <- otu_table(rowsum(x, group, reorder = reorder), taxa_are_rows = TRUE)
    stopifnot(all.equal(as.character(new_names$group), taxa_names(otu)))
    taxa_names(otu) <- new_names$archetype
    if (needs_flip) {
      otu <- t(otu)
    }
    return(otu)
  }
)

#' @rdname merge_taxa_vec
setMethod(
  "merge_taxa_vec", "taxonomyTable",
  function(x, group, reorder = FALSE, tax_adjust = 1L) {
    stopifnot(ntaxa(x) == length(group))
    # Temporary stopgap to avoid hidden errors if internal variable names are
    # in the tax table
    if (any(c(".taxon", ".group") %in% rank_names(x))) {
      stop("Currently requires that '.taxon' and '.group' are not in `rank_names(x)`")
    }
    # drop taxa with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding taxa will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    if (tax_adjust == 0L) {
      return(merge_taxa_vec_pseudo(x, group, reorder = reorder))
    } else if (tax_adjust == 1L) {
      na_bad <- FALSE
    } else if (tax_adjust == 2L) {
      na_bad <- TRUE
    }
    k <- length(rank_names(x))
    # bad_string is used to temporarily mark bad values in the tax table
    bad_string <- paste0("BAD", Sys.time())
    # Reduce each group to one row; sort if needed; then finish flushing bad
    # ranks and making new tax table
    reduced <- x %>%
      as("matrix") %>%
      as_tibble()
    reduced[, ".taxon"] <- taxa_names(x)
    reduced[, ".group"] <- factor(group, levels = unique(group))

    reduced_by_group <- as_tibble(apply(
      reduced, 2, function(xx) {
        unlist(tapply(xx, reduced$.group, bad_or_unique,
          bad = bad_string, simplify = FALSE
        ))
      }
    ))

    reduced_by_group[, ".taxon"] <-
      tapply(reduced$.taxon, reduced$.group, function(xx) {
        xx[[1]]
      })

    if (reorder) {
      reduced_by_group <- reduced_by_group %>%
        arrange(.group)
    }

    reduced_by_group <- reduced_by_group %>%
      select(-.group) %>%
      tibble::column_to_rownames(".taxon")
    # If only one tax rank, just convert bad_string -> NA; else, need to
    # propagate bad ranks downwards and convert to NAs
    if (identical(length(rank_names(x)), 1L)) {
      reduced[[1]] <- reduced[[1]] %>%
        {
          ifelse(. == bad_string, NA_character_, .)
        }
      reduced %>%
        as("matrix") %>%
        tax_table()
    } else {
      reduced %>%
        apply(1, bad_flush_right, bad = bad_string, na_bad = na_bad, k = k) %>%
        t() %>%
        tax_table()
    }
  }
)

#' @rdname merge_taxa_vec
setMethod(
  "merge_taxa_vec", "phylo",
  function(x, group) {
    merge_taxa_vec_pseudo(x, group)
  }
)

#' @rdname merge_taxa_vec
setMethod(
  "merge_taxa_vec", "XStringSet",
  function(x, group, reorder = FALSE) {
    merge_taxa_vec_pseudo(x, group, reorder = reorder)
  }
)


#' Pseudo-merge taxa in groups
#'
#' Returns `x` pruned to the first taxon of each group defined in `group`.
#'
#' @param x a phyloseq component-class object
#' @param group a vector with one element for each taxon in `x` that defines
#'   the new groups
#' @keywords internal
merge_taxa_vec_pseudo <- function(x, group, reorder = FALSE) {
  stopifnot(ntaxa(x) == length(group))
  # drop taxa with `is.na(group)`
  if (anyNA(group)) {
    warning("`group` has missing values; corresponding taxa will be dropped")
    x <- prune_taxa(!is.na(group), x)
    group <- group[!is.na(group)]
  }
  # Archetypes are the first taxon in each group
  archetypes <- tibble(
    taxon = taxa_names(x),
    group = factor(group, levels = unique(group))
  ) %>%
    group_by(group) %>%
    mutate(archetype = taxon[1])

  if (reorder) {
    archetypes %>% arrange(group)
  }
  select_taxa(x, archetypes$taxon, reorder = TRUE)
}

# helper functions ------------------------------------------------------------

#' Reduce a vector x to its unique value or the value of `bad`
#'
#' Helper for `merge_taxa_vec()`
#'
#' @param x a vector
#' @param bad the string representing a bad value
#' @keywords internal
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
bad_or_unique <- function(x, bad = "BAD") {
  if (length(unique(x)) == 1) {
    x[[1]]
  } else {
    bad
  }
}

#' Replace all values with NA upon seeing a bad value
#'
#' Helper for `merge_taxa_vec()`
#'
#' @param x a vector
#' @param bad the string representing a bad value
#' @param na_bad whether NAs should also be treated as bad
#' @param k the index to which values are flushed
#' @keywords internal
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
bad_flush_right <- function(x, bad = "BAD", na_bad = FALSE, k = length(x)) {
  if (na_bad) {
    which_bad <- which(x == bad | is.na(x))
  } else {
    which_bad <- which(x == bad)
  }
  if (length(which_bad) > 0) {
    x[seq(min(which_bad), k)] <- NA
  }
  return(x)
}

#' Merge samples by a sample variable or factor
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Firstly release in the [speedyseq](https://github.com/mikemc/speedyseq/) R
#' package by Michael R. McLaren.
#'
#' This function provides an alternative to `phyloseq::merge_samples()` that
#' better handles sample variables of different types, especially categorical
#' sample variables. It combines the samples in `x` defined by the sample
#' variable or factor `group` by summing the abundances in `otu_table(x)` and
#' combines sample variables by the summary functions in `funs`. The default
#' summary function, `unique_or_na()`, collapses the values within a group to a
#' single unique value if it exists and otherwise returns NA. The new (merged)
#' samples are named by the values in `group`.
#'
#' @param x A `phyloseq`, `otu_table`, or `sample_data` object
#' @param group A sample variable or a vector of length `nsamples(x)` defining
#'   the sample grouping. A vector must be supplied if x is an otu_table
#' @param fun_otu Function for combining abundances in the otu table; default
#'   is `sum`. Can be a formula to be converted to a function by
#'   [purrr::as_mapper()]
#' @param funs Named list of merge functions for sample variables; default is
#'   `unique_or_na`
#' @param reorder Logical specifying whether to reorder the new (merged)
#'   samples by name
#'
#' @export
#'
#' @examples
#' data(enterotype)
#'
#' # Merge samples with the same project and clinical status
#' ps <- enterotype
#' sample_data(ps) <- sample_data(ps) %>%
#'   transform(Project.ClinicalStatus = Project:ClinicalStatus)
#' sample_data(ps) %>% head()
#' ps0 <- merge_samples2(ps, "Project.ClinicalStatus",
#'   fun_otu = mean,
#'   funs = list(Age = mean)
#' )
#' sample_data(ps0) %>% head()
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X)) modified by Adrien Taudiere
setGeneric(
  "merge_samples2",
  function(x,
           group,
           fun_otu = sum,
           funs = list(),
           reorder = FALSE) {
    standardGeneric("merge_samples2")
  }
)

#' @rdname merge_samples2
setMethod(
  "merge_samples2",
  signature("phyloseq"),
  function(x, group, fun_otu = sum, funs = list(), reorder = FALSE) {
    if (length(group) == 1) {
      stopifnot(group %in% sample_variables(x))
      group <- sample_data(x)[[group]]
    } else {
      stopifnot(identical(length(group), nsamples(x)))
    }
    # Drop samples with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding samples will be dropped")
      x <- prune_samples(!is.na(group), x)
      group <- group[!is.na(group)]
    }
    # Merge
    otu.merged <- merge_samples2(otu_table(x), group,
      fun_otu = fun_otu,
      reorder = reorder
    )
    if (!is.null(access(x, "sam_data"))) {
      sam.merged <- merge_samples2(sample_data(x), group, funs = funs)
    } else {
      sam.merged <- NULL
    }
    phyloseq(
      otu.merged,
      sam.merged,
      access(x, "tax_table"),
      access(x, "phy_tree"),
      access(x, "refseq")
    )
  }
)

#' @rdname merge_samples2
setMethod(
  "merge_samples2",
  signature("otu_table"),
  function(x, group, fun_otu = sum, reorder = FALSE) {
    stopifnot(identical(length(group), nsamples(x)))
    # Work with samples as rows, and remember to flip back at end if needed
    needs_flip <- taxa_are_rows(x)
    if (needs_flip) {
      x <- t(x)
    }
    # Drop samples with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding samples will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    # Merging; result is a matrix with taxa as columns and rownames
    # corresponding to `group`
    if (identical(fun_otu, sum)) {
      x.merged <- rowsum(x, group, reorder = reorder)
    } else {
      stopifnot(!".group" %in% colnames(x))
      x.merged <- x %>%
        as("matrix") %>%
        tibble::as_tibble() %>%
        cbind(.group = group) %>%
        group_by(.group) %>%
        summarise(across(everything(), purrr::as_mapper(fun_otu)))


      if (reorder) {
        x.merged <- x.merged %>% arrange(.group)
      }
      x.merged <- x.merged %>%
        tibble::column_to_rownames(".group")
    }
    # Return an otu table in the proper orientation
    x.merged <- x.merged %>% otu_table(taxa_are_rows = FALSE)
    if (needs_flip) {
      x.merged <- t(x.merged)
    }
    return(x.merged)
  }
)

#' @rdname merge_samples2
setMethod(
  "merge_samples2",
  signature("sample_data"),
  function(x, group, funs = list(), reorder = FALSE) {
    if (length(group) == 1) {
      stopifnot(group %in% sample_variables(x))
      group <- x[[group]]
    } else {
      stopifnot(identical(length(group), nsamples(x)))
    }
    # Drop samples with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding samples will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    ## Set the functions f used to merge each sample variable.
    # Named logical vector indicating whether each variable is in the funs
    var_in_funs <- names(x) %>%
      rlang::set_names(. %in% names(funs), .)
    # For vars in the funs, run f through as_mapper; else, use the default f
    funs <- purrr::map2(
      var_in_funs, names(var_in_funs),
      ~ if (.x) purrr::as_mapper(funs[[.y]]) else unique_or_na
    )
    ## Merge variable values, creating a new sample_data object with one row
    ## per group.
    # A "sample_data" object is a list of data variables (columns); strategy is
    # to reduce each variable with `merge_groups()`, and then recombine into a
    # data.frame. The call to `merge_groups()` will sort by `group` values,
    # which we need to account for when setting the new sample names.
    new_sample_names <- group %>%
      unique() %>%
      sort() %>%
      as.character()
    x.merged <- purrr::map2(
      x, funs,
      ~ merge_groups(.x, group = group, f = .y)
    ) %>%
      data.frame() %>%
      vctrs::vec_set_names(new_sample_names)
    ## Put back in initial order
    if (!reorder) {
      initial_order <- group %>%
        unique() %>%
        as.character()
      x.merged <- x.merged[initial_order, , drop = FALSE]
    }
    ## Return as sample data with group names preserved
    x.merged %>% MiscMetabar:::sample_data_stable()
  }
)



# Helpers ---------------------------------------------------------------------

#' Get the unique value in x or NA if none
#'
#' If `unique(x)` is a single value, return it; otherwise, return an NA of the
#' same type as `x`. If `x` is a factor, then the levels and ordered status
#' will be kept in either case. If `x` is a non-atomic vector (i.e. a list),
#' then the logical `NA` will be used.
#'
#' @param x A vector
#' @export
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
#' @examples
#' f <- factor(c("a", "a", "b", "c"), ordered = TRUE)
#' unique_or_na(f)
#' unique_or_na(f[1:2])
#'
#' x <- c("a", "b", "a")
#' unique_or_na(x[c(1, 3)])
#' unique_or_na(x)
#' unique_or_na(x) %>% typeof()
unique_or_na <- function(x) {
  UseMethod("unique_or_na")
}

#' @export
unique_or_na.default <- function(x) {
  if (length(unique(x)) == 1) {
    x[[1]]
  } else if (is.atomic(x)) {
    as(NA, typeof(x))
  } else {
    NA
  }
}

#' @export
unique_or_na.factor <- function(x) {
  if (length(unique(x)) == 1) {
    x[[1]]
  } else {
    factor(NA, levels = levels(x), ordered = is.ordered(x))
  }
}

#' Merge groups of elements within a vector by a function
#'
#' Internal function used in `merge_samples2()` to merge variables. Note, owing
#' to the use of `split()`, the merged elements in the new vector will be
#' reordered according to `group`.
#'
#' @param x A vector whose elements will be merged.
#' @param group A vector such that `as.factor(group)` defines the grouping.
#' @param f A function that, when applied to a subvector of x, returns a single
#'   value. Can also be a formula as interpretted by `purrr::as_mapper()`.
#'
#' @keywords internal
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
merge_groups <- function(x, group, f = unique_or_na) {
  f <- purrr::as_mapper(f)
  split(x, group) %>%
    purrr::map(f) %>%
    {
      vctrs::vec_c(!!!., .name_spec = rlang::zap())
    }
}


#' Create sample data without adjusting row/sample names
#'
#' `phyloseq::sample_data()` will change the sample names from the row names if
#' they are `as.character(seq(1, row(object)))`. This function instead keeps the
#' names as is.
#'
#' @param object A "data.frame"-class object
#'
#' @keywords internal
#'
#' @examples
#' x <- data.frame(var1 = letters[1:3], var2 = 7:9)
#' rownames(x)
#' sample_data(x)
#' MiscMetabar:::sample_data_stable(x)
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
sample_data_stable <- function(object) {
  # Modified from phyloseq's sample_data data.frame method; see
  # https://github.com/joey711/phyloseq/blob/master/R/sampleData-class.R
  stopifnot(identical(class(object), "data.frame"))
  # Make sure there are no phantom levels in categorical variables

  object <- droplevels(as(object, "data.frame"))
  # instantiate first to check validity
  SM <- new("sample_data", object)
  return(SM)
}


# select_taxa -----------------------------------------------------------------

#' Select a subset of taxa in a specified order where possible
#'
#' Select (a subset of) taxa; if `x` allows taxa to be reordered, then taxa are
#' given in the specified order.
#'
#' This is a simple selector function that is like `prune_taxa(taxa, x)` when
#' `taxa` is a character vector but always gives the taxa in the order `taxa`
#' if possible (that is, except for phy_tree's and phyloseq objects that
#' contain phy_tree's).
#'
#' @param x A phyloseq object or phyloseq component object
#' @param taxa Character vector of taxa to select, in requested order
#' @param reorder Logical specifying whether to use the order in `taxa` (TRUE)
#'   or keep the order in `taxa_names(x)` (FALSE)
#' @keywords internal
#' @author Michael R. McLaren (orcid: [0000-0003-1575-473X](https://orcid.org/0000-0003-1575-473X))
#' @rdname select_taxa-methods
setGeneric(
  "select_taxa",
  function(x, taxa, reorder = TRUE) standardGeneric("select_taxa")
)

#' @rdname select_taxa-methods
setMethod(
  "select_taxa", signature("sample_data", "character"),
  function(x, taxa) {
    stopifnot(!anyDuplicated(taxa))
    x
  }
)

#' @rdname select_taxa-methods
setMethod(
  "select_taxa", signature("otu_table", "character"),
  function(x, taxa, reorder = TRUE) {
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder) {
      taxa <- intersect(taxa_names(x), taxa)
    }
    if (taxa_are_rows(x)) {
      x[taxa, , drop = FALSE]
    } else {
      x[, taxa, drop = FALSE]
    }
  }
)

#' @rdname select_taxa-methods
setMethod(
  "select_taxa", signature("taxonomyTable", "character"),
  function(x, taxa, reorder = TRUE) {
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder) {
      taxa <- intersect(taxa_names(x), taxa)
    }
    x[taxa, , drop = FALSE]
  }
)

#' @rdname select_taxa-methods
setMethod(
  "select_taxa", signature("XStringSet", "character"),
  function(x, taxa, reorder = TRUE) {
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder) {
      taxa <- intersect(taxa_names(x), taxa)
    }
    x[taxa]
  }
)

#' @rdname select_taxa-methods
setMethod(
  "select_taxa", signature("phylo", "character"),
  function(x, taxa) {
    # NOTE: `reorder` argument silently ignored if supplied
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    ape::keep.tip(x, taxa)
  }
)

#' @rdname select_taxa-methods
setMethod(
  "select_taxa", signature("phyloseq", "character"),
  function(x, taxa, reorder = TRUE) {
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder) {
      taxa <- intersect(taxa_names(x), taxa)
    }
    otu_table(x) <- select_taxa(otu_table(x), taxa)

    tax_order <- taxa_names(otu_table(x))
    if (!is.null(tax_table(x, FALSE))) {
      # If there is a taxonomyTable, re-order that too.
      x@tax_table <- tax_table(x)[tax_order, ]
    }
    if (!is.null(refseq(x, FALSE))) {
      # If there is a XStringSet, re-order that too.
      x@refseq <- refseq(x)[tax_order]
    }
  }
)
