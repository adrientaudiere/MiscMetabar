################################################################################
#' Resolve conflict in a vector of taxonomy values
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#'   Internally used in the function [assign_blastn()] with method="vote"
#'   and [assign_vsearch_lca()] if `top_hits_only` is FALSE and `vote_algorithm` is not NULL.
#' @param vec (required) A vector of (taxonomic) values
#' @param method One of "consensus", "rel_majority", "abs_majority",
#'   "preference" or "unanimity". See details.
#'
#' @param strict (logical, default FALSE). If TRUE, NA are considered as
#'   informative in resolving conflict (i.e. NA are taking into account in vote).
#'   See details for more informations.
#' @param second_method One of "consensus", "rel_majority", "abs_majority",
#'   or "unanimity". Only used if method = "preference". See details.
#' @param nb_agree_threshold (Int, default 1) The minimum number of times a
#'   value must arise to be selected using vote. If 2, we only kept
#'   taxonomic value present at least 2 times in the vector.
#' @param preference_index (Int. default NULL). Required if method="preference".
#'   Useless for other method. The preference index is the index of the value in
#'   vec for which we have a preference.
#' @param collapse_string (default '/'). The character to collapse taxonomic names
#'   when multiple assignment is done.
#' @param replace_collapsed_rank_by_NA (logical, default FALSE). If set to TRUE,
#'   all multiple assignments (all taxonomic rank including the 'collapse_string'
#'   parameter) are replaced by NA.
#'
#' @details
#'
#' - `unanimity`: Only keep taxonomic value when all methods are agree
#' 	 - By default, the value with NA are not taking into account (strict=FALSE)
#' 	 - If `strict` , one NA in the row is sufficient to return a NA
#'
#' - `consensus`: Keep all taxonomic values separated by a '/' (separation can be modify using param `collapse_string`)
#' 	 - If `strict` is TRUE, NA are also added to taxonomic vector such as 'Tiger/Cat/NA' instead of 'Tiger/Cat'
#'
#' - `abs_majority`: Keep the most found taxonomic value if it represent at least half of all taxonomic values.
#' 	 - If `strict` is TRUE, NA values are also count to determine the majority. For example, a vector of taxonomic rank c("A", "A", "A", "B", NA, NA) will give a value of 'A' if `strict` is FALSE (default) but a value of NA if `strict` is TRUE.
#'
#' - `rel_majority`: Keep the most found taxonomic value. In case of equality, apply a consensus strategy (collapse values separated by a '/') across the most found taxonomic values.
#' 	 - If `strict` is TRUE, NA are considered as a rank and can win the relative majority vote.  If `strict` is FALSE (default), NA are removed before ranking the taxonomic values.
#' 	 - `nb_agree_threshold`: Only keep return value when at least x methods agreed with x is set by parameter `nb_agree_threshold`. By default, (`nb_agree_threshold` = 1): a majority of one is enough.
#'
#' - `preference`: Keep the value from a preferred column.
#' 	 - when the value is NA in the preferred column, apply a second strategy (by default `consensus`) to the other column (parameter `second_method`). Note that the parameters `strict` and `nb_agree_threshold` are used for the second_method consensus.
#' @returns a vector of length 1 (one character value)
#'
#' @export
#' @author Adrien Taudière
#'
#' @examples
#'
#' resolve_vector_ranks(c("A"))
#' resolve_vector_ranks(c("A"),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A"), method = "abs_majority")
#' resolve_vector_ranks(c("A"), method = "rel_majority")
#' resolve_vector_ranks(c("A"),
#'   method = "rel_majority",
#'   nb_agree_threshold = 2
#' )
#' resolve_vector_ranks(c("A"), method = "unanimity")
#'
#' resolve_vector_ranks(c("A", "A", "A"))
#' resolve_vector_ranks(c("A", "A", "A"),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "A", "A"), method = "abs_majority")
#' resolve_vector_ranks(c("A", "A", "A"), method = "rel_majority")
#' resolve_vector_ranks(c("A", "A", "A"), method = "unanimity")
#'
#' resolve_vector_ranks(c(NA, NA, NA))
#' resolve_vector_ranks(c(NA, NA, NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c(NA, NA, NA), method = "abs_majority")
#' resolve_vector_ranks(c(NA, NA, NA), method = "rel_majority")
#' resolve_vector_ranks(c(NA, NA, NA), method = "unanimity")
#'
#' resolve_vector_ranks(c("A", "A", NA))
#' resolve_vector_ranks(c("A", "A", NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "A", NA), method = "rel_majority")
#' resolve_vector_ranks(c("A", "A", NA), method = "abs_majority")
#' resolve_vector_ranks(c("A", "A", NA, NA),
#'   method = "abs_majority",
#'   strict = FALSE
#' )
#' resolve_vector_ranks(c("A", "A", NA, NA),
#'   method = "abs_majority",
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "A", NA), method = "unanimity")
#' resolve_vector_ranks(c("A", "A", NA),
#'   method = "unanimity",
#'   strict = TRUE
#' )
#'
#' resolve_vector_ranks(c("A", "B", NA))
#' resolve_vector_ranks(c("A", "B", NA), strict = TRUE)
#' resolve_vector_ranks(c("A", "B", NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "B", NA), method = "abs_majority")
#' resolve_vector_ranks(c("A", "B", NA), method = "rel_majority")
#' resolve_vector_ranks(c("A", "B", NA),
#'   method = "rel_majority",
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "B", NA),
#'   method = "rel_majority",
#'   nb_agree_threshold = 2
#' )
#' resolve_vector_ranks(c("A", "B", NA), method = "unanimity")
#'
#' resolve_vector_ranks(c("A", NA, NA))
#' resolve_vector_ranks(c("A", NA, NA), method = "rel_majority")
#' resolve_vector_ranks(c("A", NA, NA), method = "unanimity")
#' resolve_vector_ranks(c("A", NA, NA),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", NA, NA),
#'   method = "preference",
#'   preference_index = 2
#' )
#' resolve_vector_ranks(c("A", NA, "B"),
#'   method = "preference",
#'   preference_index = 2
#' )
#' resolve_vector_ranks(c("A", NA, "B"),
#'   method = "preference",
#'   preference_index = 2, second_method = "abs_majority"
#' )
#'
#' resolve_vector_ranks(c("A", "B", "B"))
#' resolve_vector_ranks(c("A", "B", "B"),
#'   method = "preference",
#'   preference_index = 1
#' )
#' resolve_vector_ranks(c("A", "B", "B"), method = "abs_majority")
#' resolve_vector_ranks(c("A", "B", "B"), method = "rel_majority")
#' resolve_vector_ranks(c("A", "B", "B"), method = "unanimity")
#'
#'
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA))
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   method = "abs_majority"
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   method = "abs_majority",
#'   strict = TRUE
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA),
#'   method = "preference", preference_index = 6, second_method = "abs_majority"
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA, NA),
#'   method = "preference", preference_index = 6, second_method = "abs_majority"
#' )
#' resolve_vector_ranks(c("A", "A", "A", "B", NA, NA, NA),
#'   method = "preference", preference_index = 6, second_method = "abs_majority",
#'   strict = TRUE
#' )
resolve_vector_ranks <- function(vec,
                                 method = c(
                                   "consensus",
                                   "rel_majority",
                                   "abs_majority",
                                   "preference",
                                   "unanimity"
                                 ),
                                 strict = FALSE,
                                 second_method = c(
                                   "consensus",
                                   "rel_majority",
                                   "abs_majority",
                                   "unanimity"
                                 ),
                                 nb_agree_threshold = 1,
                                 preference_index = NULL,
                                 collapse_string = "/",
                                 replace_collapsed_rank_by_NA = FALSE) {
  method <- match.arg(method)
  second_method <- match.arg(second_method)

  if (sum(is.na(vec)) == length(vec)) {
    res <- NA
    return(res)
  }

  if (method == "consensus") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    res <- paste0(unique(vec), collapse = collapse_string)
  } else if (method == "preference") {
    if (is.null(preference_index)) {
      stop("You must specify a preference_index if method = 'preference'")
    }
    res <- vec[preference_index]
    if (is.na(res)) {
      res <- resolve_vector_ranks(
        vec[-preference_index],
        method = second_method,
        nb_agree_threshold = nb_agree_threshold,
        strict = strict,
        collapse_string = collapse_string,
        replace_collapsed_rank_by_NA = replace_collapsed_rank_by_NA
      )
    }
  } else if (method == "abs_majority") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    nval <- length(vec)
    higher_val <- sort(table(vec), decreasing = T)[1]
    if (higher_val / nval > 0.5) {
      res <- names(higher_val)
    } else {
      res <- NA
    }
  } else if (method == "rel_majority") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    nval <- sum(table(vec, useNA = "ifany") == max(table(vec, useNA = "ifany")))
    if (sum(sort(table(vec, useNA = "ifany"), decreasing = T)[1:nval]) >= nb_agree_threshold) {
      res <- paste0(names(sort(table(vec, useNA = "ifany"), decreasing = T)[1:nval]), collapse = collapse_string)
    } else {
      res <- NA
    }
  } else if (method == "unanimity") {
    if (!strict) {
      vec <- as.vector(na.omit(vec))
    }
    if (length(unique(vec)) == 1) {
      res <- unique(vec)
    } else {
      res <- NA
    }
  }
  if (replace_collapsed_rank_by_NA &&
    sum(grepl(collapse_string, res)) > 0) {
    res <- NA
  }

  return(res)
}
################################################################################


################################################################################
#' Format a fasta database in sintax format
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' Only tested with Unite and Eukaryome fasta file for the moment. Rely on the presence of the pattern
#'  pattern_tax default "k__" to format the header.
#'
#'  A reference database in sintax format
#'  contain taxonomic information in the header of
#'  each sequence in the form of a string starting with ";tax=" and followed
#'  by a comma-separated list of up to nine taxonomic identifiers. Each taxonomic
#'  identifier must start with an indication of the rank by one of the letters d
#'  (for domain) k (kingdom), p (phylum), c (class), o (order), f (family),
#'   g (genus), s (species), or t (strain). The letter is followed by a colon
#'    (:) and the name of that rank. Commas and semicolons are not allowed in
#'    the name of the rank. Non-ascii characters should be avoided in the names.
#'
#'  Example:
#'
#'  \>X80725_S000004313;tax=d:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Enterobacteriales,f:Enterobacteriaceae,g:Escherichia/Shigella,s:Escherichia_coli,t:str._K-12_substr._MG1655
#'
#' @param fasta_db A link to a fasta files
#' @param taxnames A list of names to format. You must specify either fasta_db OR taxnames, not both.
#' @param pattern_tax (default "k__") The pattern to replace by pattern_sintax.
#' @param pattern_sintax (default "tax=k:") Useless for most users. Sometimes you may want to
#'   replacte by "tax=d:" (d for domain instead of kingdom).
#' @param output_path (optional) A path to an output fasta files. Only used if fasta_db is set.
#' @export
#' @author Adrien Taudière
#' @seealso [format2dada2_species()], [format2dada2()]
#' @return Either an object of class DNAStringSet or a vector of reformated names
format2sintax <- function(fasta_db = NULL,
                          taxnames = NULL,
                          pattern_tax = "k__",
                          pattern_sintax = "tax=k:",
                          output_path = NULL) {
  if (is.null(taxnames) && is.null(fasta_db)) {
    stop("You must specify taxnames or fasta_db parameter.")
  } else if (!is.null(taxnames) && !is.null(fasta_db)) {
    stop("You must specify either taxnames or fasta_db, not both.")
  } else if (!is.null(taxnames)) {
    new_names <- taxnames %>%
      {
        gsub(";", ",", .)
      } %>%
      {
        gsub(pattern_k, paste0(";", pattern_sintax), .)
      } %>%
      {
        gsub("__", ":", .)
      } %>%
      {
        gsub(";;", ";", .)
      } %>%
      {
        gsub(paste0(",", pattern_sintax), paste0(";", pattern_sintax), .)
      }
    return(new_names)
  } else if (!is.null(fasta_db)) {
    dna <- Biostrings::readDNAStringSet(fasta_db)
    new_names <- names(dna) %>%
      {
        gsub(";", ",", .)
      } %>%
      {
        gsub(pattern_tax, paste0(";", pattern_sintax), .)
      } %>%
      {
        gsub("__", ":", .)
      } %>%
      {
        gsub(";;", ";", .)
      } %>%
      {
        gsub(paste0(",", pattern_sintax), paste0(";", pattern_sintax), .)
      }

    names(dna) <- new_names
    if (!is.null(output_path)) {
      Biostrings::writeXStringSet(dna, filepath = output_path)
    }
    return(dna)
  }
}
################################################################################


################################################################################
#' Format a fasta database in dada2 format
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' First format in sintax format and then in dada2 format
#'
#' @param fasta_db A link to a fasta files
#' @param taxnames A list of names to format. You must specify either fasta_db OR taxnames, not both.
#' @param output_path (optional) A path to an output fasta files. Only used if fasta_db is set.
#' @param from_sintax (logical, default FALSE) Is the original fasta file in sintax format?
#' @param pattern_to_remove (a regular expression) Define a pattern to remove. For example,
#'   pattern_to_remove = "\\|rep.*" remove all character after '|rep' to force [dada2::assignTaxonomy()]
#'   to not use the database as a Unite-formated database
#' @param ... Additional arguments passed on to [format2sintax()] function
#' @export
#' @author Adrien Taudière
#' @return Either an object of class DNAStringSet or a vector of reformated names
#' @seealso [format2dada2_species()], [format2sintax()]

format2dada2 <- function(fasta_db = NULL,
                         taxnames = NULL,
                         output_path = NULL,
                         from_sintax = TRUE,
                         pattern_to_remove = NULL,
                         ...) {
  if (is.null(taxnames) && is.null(fasta_db)) {
    stop("You must specify taxnames or fasta_db parameter.")
  } else if (!is.null(taxnames) && !is.null(fasta_db)) {
    stop("You must specify either names or fasta_db, not both.")
  } else if (!is.null(taxnames)) {
    if (from_sintax) {
      new_names <- taxnames
    } else {
      new_names <- format2sintax(taxnames = taxnames, ...)
    }
    new_names <- new_names |>
      stringr::str_split_fixed(";tax=", n = 2) |>
      as_tibble() |>
      tidyr::unite(taxnames, c(V2, V1), sep = ";") |>
      pull(taxnames) %>%
      {
        gsub(":", "__", .)
      } %>%
      {
        gsub(",", ";", .)
      }

    if (!is.null(pattern_to_remove)) {
      new_names <- new_names |>
        stringr::str_remove(pattern_to_remove)
    }

    return(new_names)
  } else if (!is.null(fasta_db)) {
    dna <- Biostrings::readDNAStringSet(fasta_db)
    if (from_sintax) {
      new_names <- names(dna)
    } else {
      new_names <- format2sintax(names(dna), ...)
    }

    # Add the good number of level to each line
    new_names <- map_chr(new_names, ~ {
      nb_char <- str_count(.x, ":")
      diff <- max_char - nb_char - 2
      if (diff > 0) {
        return(paste0(.x, strrep(",", diff)))
      } else {
        return(.x)
      }
    })


    new_names <- new_names |>
      stringr::str_split_fixed(";tax=", n = 2) |>
      as_tibble() |>
      tidyr::unite(taxnames, c(V2, V1), sep = "") |>
      pull(taxnames) |>
      paste0(";") |>
      gsub(pattern = ":", replacement = "__") |>
      gsub(pattern = ",", replacement = ";")


    if (!is.null(pattern_to_remove)) {
      new_names <- new_names |>
        stringr::str_remove(pattern_to_remove)
    }

    names(dna) <- new_names

    if (!is.null(output_path)) {
      Biostrings::writeXStringSet(dna, filepath = output_path)
      invisible(dna)
    } else {
      return(dna)
    }
  }
}
################################################################################


################################################################################
#' Format a fasta database in dada2 format for Species assignment
#'
#' @description
#'
#' <a href="https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle">
#' <img src="https://img.shields.io/badge/lifecycle-experimental-orange" alt="lifecycle-experimental"></a>
#'
#' First format in sintax format and then in dada2 format
#'
#' @param fasta_db A link to a fasta files
#' @param taxnames A list of names to format. You must specify either fasta_db OR taxnames, not both.
#' @param output_path (optional) A path to an output fasta files. Only used if fasta_db is set.
#' @param from_sintax (logical, default FALSE) Is the original fasta file in sintax format?
#' @param ... Additional arguments passed on to [format2sintax()] function
#' @export
#' @author Adrien Taudière
#' @return Either an object of class DNAStringSet or a vector of reformated names
#' @seealso [format2dada2_species()], [format2sintax()]
#'
format2dada2_species <- function(
  fasta_db = NULL,
  taxnames = NULL,
  from_sintax = FALSE,
  output_path = NULL,
  ...
) {
  if (is.null(taxnames) && is.null(fasta_db)) {
    stop("You must specify taxnames or fasta_db parameter.")
  } else if (!is.null(taxnames) && !is.null(fasta_db)) {
    stop("You must specify either taxnames or fasta_db, not both.")
  } else if (!is.null(taxnames)) {
    if (from_sintax) {
      new_names <- paste(
        stringr::str_extract(taxnames, "^(.*?);tax=", group = T),
        stringr::str_extract(taxnames, "g:(.*?),", group = T),
        stringr::str_extract(taxnames, "s:(.*?)$", group = T),
        sep = " "
      )
    } else {
      new_names <- paste(
        stringr::str_extract(taxnames, "^(.*?)k__", group = T),
        stringr::str_extract(taxnames, "g__(.*?);", group = T),
        stringr::str_extract(taxnames, "s__(.*?)$", group = T),
        sep = " "
      )
    }
    return(new_names)
  } else if (!is.null(fasta_db)) {
    dna <- Biostrings::readDNAStringSet(fasta_db)
    taxnames <- names(dna)
    if (from_sintax) {
      id <- stringr::str_extract(taxnames, "^(.*?);tax=", group = T)
      id[is.na(id)] <- taxnames[is.na(id)]
      genus <- stringr::str_extract(taxnames, "g:(.*?),", group = T)
      species <- stringr::str_extract(taxnames, "s:(.*?)$", group = T)
    } else {
      id <- stringr::str_extract(taxnames, "^(.*?)k__", group = T)
      id[is.na(id)] <- taxnames[is.na(id)]
      genus <- stringr::str_extract(taxnames, "g__(.*?);", group = T)
      species <- stringr::str_extract(taxnames, "s__(.*?)$", group = T)
    }
    new_names <- paste(id, genus, species, sep = " ")

    names(dna) <- new_names
    names(dna)[is.na(names(dna))] <- taxnames[names(dna)]

    if (!is.null(output_path)) {
      Biostrings::writeXStringSet(dna, filepath = output_path)
    }
    return(dna)
  }
}
################################################################################
