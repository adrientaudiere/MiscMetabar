################################################################################
#' Make a datatable with the taxonomy of a \code{\link{phyloseq-class}} object
#' @description
#' `r lifecycle::badge("maturing")`
#' @inheritParams clean_pq
#' @param abundance (default: TRUE) Does the number of sequences is print
#' @param taxonomic_level (default: NULL) a vector of selected taxonomic
#' level using their column numbers (e.g. taxonomic_level = c(1:7))
#' @param modality (default: NULL) A sample modality to split
#' OTU abundancy by level of the modality
#' @param ... Other argument for the datatable function
#'
#'
#' @author Adrien Taudière
#' @return A datatable
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' tax_datatable(subset_taxa(
#'   GlobalPatterns,
#'   rowSums(GlobalPatterns@otu_table) > 10000
#' ))
#'
#' # Using modality
#' tax_datatable(GlobalPatterns,
#'   modality = GlobalPatterns@sam_data$SampleType
#' )
tax_datatable <- function(physeq,
                          abundance = TRUE,
                          taxonomic_level = NULL,
                          modality = NULL,
                          ...) {
  df <- as.data.frame(unclass(physeq@tax_table))

  if (!is.null(taxonomic_level)) {
    df <- df[, taxonomic_level]
  }

  if (is.null(modality)) {
    if (abundance) {
      if (physeq@otu_table@taxa_are_rows) {
        df$nb_seq <- rowSums(physeq@otu_table)
      } else {
        df$nb_seq <- colSums(physeq@otu_table)
      }
    }
  } else {
    modality <- as.factor(modality)
    if (physeq@otu_table@taxa_are_rows) {
      for (mod in levels(modality)) {
        varname <- paste0("nb_seq_", mod)
        df[[varname]] <- rowSums(physeq@otu_table[, modality == mod])
      }
      df$nb_seq_tot <- rowSums(physeq@otu_table)
    } else {
      for (mod in levels(modality)) {
        varname <- paste0("nb_seq_", mod)
        df[[varname]] <- colSums(physeq@otu_table[modality == mod, ])
      }
      df$nb_seq_tot <- colSums(physeq@otu_table)
    }
  }

  if (is.null(modality)) {
    dt <- DT::datatable(df, ...) %>% DT::formatStyle(
      "nb_seq",
      background = DT::styleColorBar(df$nb_seq, "steelblue"),
      backgroundSize = "100% 90%",
      backgroundRepeat = "no-repeat",
      backgroundPosition = "center"
    )
  } else {
    dt <- DT::datatable(df, ...)
    for (cn in colnames(df)[grepl("nb_seq", colnames(df))]) {
      dt <- dt %>%
        DT::formatStyle(
          cn,
          background = DT::styleColorBar(df[[cn]], "steelblue"),
          backgroundSize = "90% 90%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        )
    }
  }
  return(dt)
}
################################################################################

################################################################################
#' Compare samples in pairs using diversity and number of ASV including
#' shared ASV.
#' @description
#' `r lifecycle::badge("experimental")` #'   For the moment refseq slot need to be not Null.
#'
#' @inheritParams clean_pq
#' @param bifactor (required) a factor (present in the `sam_data` slot of
#'   the physeq object) presenting the pair names
#' @param modality the name of the column in the `sam_data`
#'   slot of the physeq object to split samples by pairs
#' @param merge_sample_by a vector to determine
#'   which samples to merge using the
#'   \code{\link[speedyseq]{merge_samples2}} function.
#'   Need to be in \code{physeq@sam_data}
#' @param nb_min_seq minimum number of sequences per sample
#'   to count the ASV/OTU
#' @param vegIndex (default: "shannon") index for the `vegan::diversity` function
#' @return A tibble
#' @importFrom rlang .data
#' @export
#' @examples
#' data(data_fungi)
#' data_fungi_low_high <- subset_samples(data_fungi, "Height" %in% c("Low", "High"))
#' compare_pairs_pq(data_fungi_low_high, bifactor = "Height", merge_sample_by = "Height")
compare_pairs_pq <- function(physeq = NULL,
                             bifactor = NULL,
                             modality = NULL,
                             merge_sample_by = NULL,
                             nb_min_seq = 0,
                             vegIndex = "shannon") {
  physeq <- clean_pq(physeq,
    clean_samples_names = FALSE,
    force_taxa_as_columns = TRUE,
    silent = T
  )

  if (!is.null(merge_sample_by)) {
    physeq <- speedyseq::merge_samples2(physeq, merge_sample_by)
    physeq <- clean_pq(physeq)
  }

  if (!is.factor(physeq@sam_data[[bifactor]])) {
    physeq@sam_data[[bifactor]] <- as.factor(physeq@sam_data[[bifactor]])
  }
  if (nlevels(as.factor(physeq@sam_data[[bifactor]]) != 2)) {
    stop("The bifactor arguments needs only two levels")
  }
  lev1 <- levels(physeq@sam_data[[bifactor]])[1]
  lev2 <- levels(physeq@sam_data[[bifactor]])[2]

  res <- list()
  if (!is.null(modality)) {
    physeq@sam_data[[modality]] <- as.factor(physeq@sam_data[[modality]])
    nmodality <- levels(physeq@sam_data[[modality]])
  } else {
    nmodality <- 1
  }
  # else {
  #   physeq@sam_data[[modality]] <- rep("mod", nsamples(physeq))
  # }

  for (i in nmodality) {
    newphyseq <- physeq
    if (!is.null(modality)) {
      newDF <- newphyseq@sam_data[newphyseq@sam_data[[modality]] == i, ]
      sample_data(newphyseq) <- sample_data(newDF)
    }
    if (nsamples(newphyseq) != 2) {
      res[[i]] <- c(NA, NA, NA, NA, NA)
      warning("At least one case do not contain 2 samples, so NA were introduced.")
    } else {
      cond1 <- newphyseq@sam_data[[bifactor]] == lev1
      cond2 <- newphyseq@sam_data[[bifactor]] == lev2
      nb_first <- rowSums(newphyseq@otu_table[cond1, ] > nb_min_seq)
      nb_second <- rowSums(newphyseq@otu_table[cond2, ] > nb_min_seq)
      nb_shared <- rowSums(newphyseq@otu_table[cond1, ] > nb_min_seq &
        newphyseq@otu_table[cond2, ] > nb_min_seq)

      div_first <- round(vegan::diversity(newphyseq@otu_table,
        index = vegIndex
      )[cond1], 2)
      div_second <- round(vegan::diversity(newphyseq@otu_table,
        index = vegIndex
      )[cond2], 2)

      res[[i]] <- c(nb_first, nb_second, nb_shared, div_first, div_second)
    }
  }

  res_df <- t(data.frame(res))

  res_df <- tibble::as_tibble(res_df)
  res_df <- res_df %>%
    mutate(percent_shared_lv1 = round(100 * .data$V3 /
      .data$V1, 2)) %>%
    mutate(percent_shared_lv2 = round(100 * .data$V3 /
      .data$V2, 2)) %>%
    mutate(ratio_nb_lv1_lv2 = round(.data$V1 /
      .data$V2, 3)) %>%
    mutate(ratio_div_lv1_lv2 = round(.data$V4 /
      .data$V5, 3))

  colnames(res_df) <- c(
    paste0("nb_", lev1),
    paste0("nb_", lev2),
    "shared",
    paste0("div_", lev1),
    paste0("div_", lev2),
    paste0("percent_shared_", lev1),
    paste0("percent_shared_", lev2),
    paste0("ratio_nb_", lev1, "_", lev2),
    paste0("ratio_div_", lev1, "_", lev2)
  )

  res_df[[modality]] <- names(res)
  res_df <- res_df %>%
    dplyr::filter(!is.na(nb_shared)) %>%
    relocate(all_of(modality))

  return(res_df)
}
