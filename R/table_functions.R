################################################################################
#' Make a datatable with the taxonomy of a \code{\link{phyloseq-class}} object
#' @description
<a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"><img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#' @inheritParams clean_pq
#' @param abundance (default: TRUE) Does the number of sequences is print
#' @param taxonomic_level (default: NULL) a vector of selected taxonomic
#' level using their column numbers (e.g. taxonomic_level = 1:7)
#' @param modality (default: NULL) A sample modality to split
#' OTU abundancy by level of the modality
#' @param ... Other argument for the datatable function
#'
#'
#' @author Adrien Taudière
#' @return A datatable
#' @export
#'
#' @examplesIf tolower(Sys.info()[["sysname"]]) != "windows"
#' data("GlobalPatterns", package = "phyloseq")
#' if (requireNamespace("DT")) {
#'   tax_datatable(subset_taxa(
#'     GlobalPatterns,
#'     rowSums(GlobalPatterns@otu_table) > 10000
#'   ))
#'
#'   # Using modality
#'   tax_datatable(GlobalPatterns,
#'     modality = GlobalPatterns@sam_data$SampleType
#'   )
#' }
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
<a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"><img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a> #'   For the moment refseq slot need to be not Null.
#'
#' @inheritParams clean_pq
#' @param bifactor (required) a factor (present in the `sam_data` slot of
#'   the physeq object) presenting the pair names
#' @param modality the name of the column in the `sam_data`
#'   slot of the physeq object to split samples by pairs
#' @param merge_sample_by a vector to determine
#'   which samples to merge using the
#'   [merge_samples2()] function.
#'   Need to be in \code{physeq@sam_data}
#' @param nb_min_seq minimum number of sequences per sample
#'   to count the ASV/OTU
#' @param veg_index (default: "shannon") index for the `vegan::diversity` function
#' @param na_remove (logical, default TRUE) If set to TRUE, remove samples with
#'   NA in the variables set in bifactor, modality and merge_sample_by.
#'   NA in variables are well managed even if na_remove = FALSE, so na_remove may
#'   be useless.
#' @return A tibble with information about the number of shared ASV, shared number of sequences
#'   and diversity
#' @importFrom rlang .data
#' @export
#' @examples
#' data_fungi_low_high <- subset_samples(data_fungi, Height %in% c("Low", "High"))
#' compare_pairs_pq(data_fungi_low_high, bifactor = "Height", merge_sample_by = "Height")
#' compare_pairs_pq(data_fungi_low_high,
#'   bifactor = "Height",
#'   merge_sample_by = "Height", modality = "Time"
#' )
compare_pairs_pq <- function(physeq = NULL,
                             bifactor = NULL,
                             modality = NULL,
                             merge_sample_by = NULL,
                             nb_min_seq = 0,
                             veg_index = "shannon",
                             na_remove = TRUE) {
  physeq <- taxa_as_columns(physeq)

  if (na_remove) {
    new_physeq <- subset_samples_pq(physeq, !is.na(physeq@sam_data[[bifactor]]))
    if (nsamples(physeq) - nsamples(new_physeq) > 0) {
      message(paste0(
        nsamples(physeq) - nsamples(new_physeq),
        " were discarded due to NA in variable bifactor."
      ))
    }
    physeq <- new_physeq
    if (!is.null(merge_sample_by)) {
      new_physeq <- subset_samples_pq(physeq, !is.na(physeq@sam_data[[merge_sample_by]]))
      if (nsamples(physeq) - nsamples(new_physeq) > 0) {
        message(paste0(
          nsamples(physeq) - nsamples(new_physeq),
          " were discarded due to NA in variable merge_sample_by."
        ))
      }
      physeq <- new_physeq
    }

    if (!is.null(modality)) {
      new_physeq <- subset_samples_pq(physeq, !is.na(physeq@sam_data[[modality]]))
      if (nsamples(physeq) - nsamples(new_physeq) > 0) {
        message(paste0(
          nsamples(physeq) - nsamples(new_physeq),
          " samples were discarded due to NA in variable modality."
        ))
      }
      physeq <- new_physeq
    }
  }

  if (!is.null(merge_sample_by)) {
    if (is.null(modality)) {
      physeq <- merge_samples2(physeq, merge_sample_by)
    } else {
      physeq@sam_data[["merge_sample_by___modality"]] <- paste0(physeq@sam_data[[merge_sample_by]], " - ", physeq@sam_data[[modality]])
      physeq <- merge_samples2(physeq, "merge_sample_by___modality")
    }
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
    nmodality <- bifactor
  }

  for (i in nmodality) {
    newphyseq <- physeq
    if (!is.null(modality)) {
      new_DF <- newphyseq@sam_data[newphyseq@sam_data[[modality]] == i, ]
      sample_data(newphyseq) <- sample_data(new_DF)
    }
    if (nsamples(newphyseq) != 2) {
      res[[i]] <- rep(NA, 8)
      warning("At least one case do not contain 2 samples, so NA were introduced.")
    } else {
      cond1 <- newphyseq@sam_data[[bifactor]] == lev1
      cond2 <- newphyseq@sam_data[[bifactor]] == lev2
      nb_first <- rowSums(newphyseq@otu_table[cond1, ] > nb_min_seq)
      nb_second <- rowSums(newphyseq@otu_table[cond2, ] > nb_min_seq)
      nb_shared <- rowSums(newphyseq@otu_table[cond1, ] > nb_min_seq &
        newphyseq@otu_table[cond2, ] > nb_min_seq)

      div_first <- round(vegan::diversity(newphyseq@otu_table,
        index = veg_index
      )[cond1], 2)
      div_second <- round(vegan::diversity(newphyseq@otu_table,
        index = veg_index
      )[cond2], 2)

      nb_shared_seq <- sum(newphyseq@otu_table[, newphyseq@otu_table[cond1, ] > nb_min_seq &
        newphyseq@otu_table[cond2, ] > nb_min_seq])

      perc_seq_shared_lv1 <- round(100 * nb_shared_seq / sum(newphyseq@otu_table[, newphyseq@otu_table[cond1, ] > nb_min_seq]), 2)

      perc_seq_shared_lv2 <- round(100 * nb_shared_seq / sum(newphyseq@otu_table[, newphyseq@otu_table[cond2, ] > nb_min_seq]), 2)

      res[[i]] <- c(nb_first, nb_second, nb_shared, div_first, div_second, nb_shared_seq, perc_seq_shared_lv1, perc_seq_shared_lv2)
    }
  }

  res_df_t <- t(as_tibble(res, .name_repair = "minimal"))
  colnames(res_df_t) <- paste0("V", seq_len(ncol(res_df_t)))
  res_df <- as_tibble(res_df_t)
  # res_df <- as_tibble(t(as_tibble(res, .name_repair = "universal")), .name_repair = "universal")

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
    paste0("nb_ASV_", lev1),
    paste0("nb_ASV_", lev2),
    "nb_shared_ASV",
    paste0("div_", lev1),
    paste0("div_", lev2),
    "nb_shared_seq",
    paste0("percent_shared_seq_", lev1),
    paste0("percent_shared_seq_", lev2),
    paste0("percent_shared_ASV_", lev1),
    paste0("percent_shared_ASV_", lev2),
    paste0("ratio_nb_", lev1, "_", lev2),
    paste0("ratio_div_", lev1, "_", lev2)
  )

  res_df$modality <- names(res)
  res_df <- res_df %>%
    dplyr::filter(!is.na(nb_shared)) %>%
    relocate(modality)

  return(res_df)
}
################################################################################


################################################################################
#' Create an visualization table to describe taxa distribution across a modality
#' @description
<a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle"><img src="https://img.shields.io/badge/lifecycle-maturing-blue" alt="lifecycle-maturing"></a>
#' @inheritParams clean_pq
#' @param modality (required) The name of a column present in the `@sam_data` slot
#'   of the physeq object. Must be a character vector or a factor.
#' @param taxonomic_levels (default = c("Phylum", "Order", "Family", "Genus"))
#'   The taxonomic levels (must be present in the `@sam_data` slot) you want to
#'   see and/or used (for example to compute a color) in the table.
#' @param min_nb_seq_taxa (default = 1000) filter out taxa with less than `min_nb_seq_taxa`
#'    sequences
#' @param log10trans (logical, default TRUE) Do sequences count is log10 transformed
#'   (using log10(x + 1) to allow 0)
#' @param void_style (logical, default FALSE) Do the default style is discard ?
#' @param lev_col_taxa Taxonomic level used to plot the background color of taxa names
#' @param arrange_by The column used to sort the table.
#'   Can take the values NULL, "proportion_samp", "nb_seq" (default), , "nb_sam" "OTU", or a column names
#'   from the levels of modality or from taxonomic levels
#' @param descending_order (logical, default TRUE) Do we use descending order when sort the table
#'  (if arrange_by is not NULL) ?
#' @param na_remove (logical, default TRUE) if TRUE remove all the samples
#'   with NA in the `split_by` variable of the `physeq@sam_data` slot
#' @param formattable_args Other args to the formattable function. See examples and `formattable::formattable()`
#'
#' @seealso `formattable::formattable()`
#'
#' @author Adrien Taudière
#' @return A datatable
#' @export
#' @importFrom grDevices col2rgb
#' @importFrom stats runif
#' @examples
#' if (requireNamespace("formattable")) {
#'   ## Distribution of the nb of sequences per OTU across Height
#'   ##   modality (nb of sequences are log-transformed).
#'   ## Only OTU with more than 10000 sequences are taking into account
#'   ## The Phylum column is discarded
#'   formattable_pq(
#'     data_fungi,
#'     "Height",
#'     min_nb_seq_taxa = 10000,
#'     formattable_args = list("Phylum" = FALSE),
#'     log10trans = TRUE
#'   )
#'
#'   ## Distribution of the nb of samples per OTU across Height modality
#'   ## Only OTU  present in more than 50 samples are taking into account
#'   formattable_pq(
#'     as_binary_otu_table(data_fungi),
#'     "Height",
#'     min_nb_seq_taxa = 50,
#'     formattable_args = list("nb_seq" = FALSE),
#'   )
#'
#'   ## Distribution of the nb of sequences per OTU across Time modality
#'   ##  arranged by Family Name in ascending order.
#'   ## Only OTU with more than 10000 sequences are taking into account
#'   ## The Phylum column is discarded
#'   formattable_pq(
#'     data_fungi,
#'     "Time",
#'     min_nb_seq_taxa = 10000,
#'     taxonomic_levels = c("Order", "Family", "Genus", "Species"),
#'     formattable_args = list(
#'       Order = FALSE,
#'       Species = formattable::formatter(
#'         "span",
#'         style = x ~ formattable::style(
#'           "font-style" = "italic",
#'           `color` = ifelse(is.na(x), "white", "grey")
#'         )
#'       )
#'     ),
#'     arrange_by = "Family",
#'     descending_order = FALSE
#'   )
#' }
#' \donttest{
#' if (requireNamespace("formattable")) {
#'   ## Distribution of the nb of sequences per OTU across Height modality
#'   ##  (nb of sequences are log-transformed).
#'   ## OTU name background is light gray for Basidiomycota
#'   ##  and dark grey otherwise (Ascomycota)
#'   ## A different color is defined for each modality level
#'   formattable_pq(
#'     data_fungi,
#'     "Height",
#'     taxonomic_levels = c("Phylum", "Family", "Genus"),
#'     void_style = TRUE,
#'     formattable_args = list(
#'       OTU = formattable::formatter(
#'         "span",
#'         style = ~ formattable::style(
#'           "display" = "block",
#'           `border-radius` = "5px",
#'           `background-color` = ifelse(Phylum == "Basidiomycota", transp("gray"), "gray")
#'         ),
#'         `padding-right` = "2px"
#'       ),
#'       High = formattable::formatter(
#'         "span",
#'         style = x ~ formattable::style(
#'           "font-size" = "80%",
#'           "display" = "inline-block",
#'           direction = "rtl",
#'           `border-radius` = "0px",
#'           `padding-right` = "2px",
#'           `background-color` = formattable::csscolor(formattable::gradient(
#'             as.numeric(x), transp("#1a91ff"), "#1a91ff"
#'           )),
#'           width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
#'         )
#'       ),
#'       Low = formattable::formatter(
#'         "span",
#'         style = x ~ formattable::style(
#'           "font-size" = "80%",
#'           "display" = "inline-block",
#'           direction = "rtl",
#'           `border-radius` = "0px",
#'           `padding-right` = "2px",
#'           `background-color` = formattable::csscolor(formattable::gradient(
#'             as.numeric(x),
#'             transp("green"), "green"
#'           )),
#'           width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
#'         )
#'       ),
#'       Middle = formattable::formatter(
#'         "span",
#'         style = x ~ formattable::style(
#'           "font-size" = "80%",
#'           "display" = "inline-block",
#'           direction = "rtl",
#'           `border-radius` = "0px",
#'           `padding-right` = "2px",
#'           `background-color` = formattable::csscolor(formattable::gradient(
#'             as.numeric(x), transp("orange"), "orange"
#'           )),
#'           width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
#'         )
#'       )
#'     )
#'   )
#' }
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `formattable::formattable()` if you
#'   use this function.
formattable_pq <- function(physeq,
                           modality,
                           taxonomic_levels = c("Phylum", "Order", "Family", "Genus"),
                           min_nb_seq_taxa = 1000,
                           log10trans = FALSE,
                           void_style = FALSE,
                           lev_col_taxa = "Phylum",
                           arrange_by = "nb_seq",
                           descending_order = TRUE,
                           na_remove = TRUE,
                           formattable_args = NULL) {
  verify_pq(physeq)

  if (na_remove) {
    new_physeq <-
      subset_samples_pq(physeq, !is.na(physeq@sam_data[[modality]]))
    if (nsamples(physeq) - nsamples(new_physeq) > 0) {
      message(paste0(
        nsamples(physeq) - nsamples(new_physeq),
        " samples were discarded due to NA in variable modality"
      ))
    }
  }
  new_physeq2 <-
    clean_pq(subset_taxa_pq(new_physeq, taxa_sums(new_physeq) > min_nb_seq_taxa))

  psm <- psmelt(new_physeq2)

  if (log10trans) {
    psm2 <- psm %>%
      group_by_at(c(modality, "OTU", taxonomic_levels)) %>%
      summarise(Ab = round(log10(1 + sum(Abundance)), 2)) %>%
      tidyr::spread(modality, Ab)

    psm3 <-
      left_join(
        psm2,
        data.frame(
          "OTU" = taxa_names(physeq),
          "proportion_samp" = round(
            taxa_sums(as_binary_otu_table(physeq)) / nsamples(physeq),
            2
          ),
          "nb_seq" = round(log10(taxa_sums(physeq)), 2),
          "nb_sam" = round(taxa_sums(as_binary_otu_table(physeq)), 2)
        )
      )
  } else {
    psm2 <- psm %>%
      group_by_at(c(modality, "OTU", taxonomic_levels)) %>%
      summarise(Ab = sum(Abundance)) %>%
      tidyr::spread(modality, Ab)

    psm3 <-
      left_join(
        psm2,
        data.frame(
          "OTU" = taxa_names(physeq),
          "proportion_samp" = round(
            taxa_sums(as_binary_otu_table(physeq)) / nsamples(physeq),
            2
          ),
          "nb_seq" = taxa_sums(physeq),
          "nb_sam" = round(taxa_sums(as_binary_otu_table(physeq)), 2)
        )
      )
  }
  if (!is.null(arrange_by)) {
    if (descending_order) {
      psm3 <- psm3 %>% arrange(desc(.data[[arrange_by]]))
    } else {
      psm3 <- psm3 %>% arrange(.data[[arrange_by]])
    }
  }
  if (void_style) {
    ftab <- formattable::formattable(psm3, formattable_args)
    return(ftab)
  } else {
    color_back_OTU <- fac2col(psm3[[lev_col_taxa]], na.col = "grey")
    formattable_args <- c(
      formattable_args,
      list(
        OTU = formattable::formatter(
          "span",
          style = ~ formattable::style(
            "display" = "block",
            `border-radius` = "5px",
            `background-color` = color_back_OTU
          ),
          `padding-right` = "2px"
        ),
        formattable::area(col = levels(as.factor(
          new_physeq2@sam_data[[modality]]
        ))) ~ formattable::formatter(
          "span",
          style = x ~ formattable::style(
            "font-size" = "80%",
            "display" = "inline-block",
            direction = "rtl",
            `border-radius` = "0px",
            `padding-right` = "2px",
            `background-color` = ifelse(x == 0, "white", formattable::csscolor(
              formattable::gradient(as.numeric(x), transp("#1a9641ff"), "#1a9641ff")
            )),
            width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
          )
        ),
        Family = formattable::formatter(
          "span",
          style = ~ formattable::style(
            "display" = "block",
            `border-radius` = "5px",
            `background-color` = formattable::csscolor(transp(fac2col(Family)))
          ),
          `padding-right` = "2px"
        ),
        Order = formattable::formatter(
          "span",
          style = ~ formattable::style(
            "display" = "block",
            `border-radius` = "5px",
            `background-color` = formattable::csscolor(transp(fac2col(Order, col.pal = viridis::viridis_pal()), 0.7))
          ),
          `padding-right` = "2px"
        ),
        Genus = formattable::formatter("span", style = x ~ formattable::style("font-style" = "italic")),
        nb_seq = formattable::formatter(
          "span",
          style = x ~ formattable::style(
            "font-size" = "80%",
            "display" = "inline-block",
            direction = "rtl",
            `border-radius` = "0px",
            `padding-right` = "5px",
            `background-color` = ifelse(x == 0, "white", formattable::csscolor(
              formattable::gradient(as.numeric(x), transp("#4d4888ff"), "#4d4888ff")
            )),
            width = formattable::percent(formattable::proportion(as.numeric(x), na.rm = TRUE))
          )
        ),
        proportion_samp = formattable::formatter(
          "span",
          style = x ~ formattable::style(
            "font-size" = "80%",
            "display" = "inline-block",
            direction = "rtl",
            `border-radius` = "0px",
            `padding-right` = "5px",
            `background-color` = ifelse(x == 0, "white", formattable::csscolor(
              formattable::gradient(as.numeric(x), transp("#1f78b4ff"), "#1f78b4ff")
            )),
            width = formattable::percent(as.numeric(x))
          )
        ),
        nb_sam = FALSE
      )
    )
    ftab <- formattable::formattable(psm3, formattable_args)
    return(ftab)
  }
}
################################################################################
