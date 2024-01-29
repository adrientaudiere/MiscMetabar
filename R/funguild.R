#' Retrieve the FUNGuild database
#' @description
#' `r lifecycle::badge("stable")`
#' The original function and documentation was written by Brendan Furneaux
#' in the [FUNGuildR](https://github.com/brendanf/FUNGuildR/) package.
#'
#' Please cite this [publication](https://www.sciencedirect.com/science/article/abs/pii/S1754504815000847).
#'
#' @param db_url a length 1 character string giving the URL to retrieve the database
#'     from
#'
#' @return a [`tibble::tibble`] containing the database, which can be passed
#'     to the `db` argument of [funguild_assign()]
#' @export
#'
#' @examples
#' get_funguild_db()
#' @references Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J,
#' Schilling JS, Kennedy PG. 2016. *FUNGuild: An open annotation tool for
#' parsing fungal community datasets by ecological guild*. Fungal Ecology
#' 20:241-248.
#' @author Brendan Furneaux (orcid: [0000-0003-3522-7363](https://orcid.org/0000-0003-3522-7363)),
#' modified by Adrien Taudière
#'
get_funguild_db <- function(db_url = "http://www.stbates.org/funguild_db_2.php") {
  taxon <- NULL
  httr::GET(url = db_url) %>%
    httr::content(as = "text") %>%
    stringr::str_split("\n") %>%
    unlist() %>%
    magrittr::extract(7) %>%
    stringr::str_replace("^\\[", "") %>%
    stringr::str_replace("]</body>$", "") %>%
    stringr::str_replace_all("\\} ?, ?\\{", "} \n {") %>%
    stringr::str_split("\n") %>%
    unlist() %>%
    purrr::map_dfr(
      function(record) {
        current_record <- jsonlite::fromJSON(record)
        if (!is.null(current_record[["TrophicMode"]])) {
          current_record$trophicMode <- current_record$TrophicMode
        }
        if (!is.null(current_record[["growthMorphology"]])) {
          current_record$growthForm <- current_record$growthMorphology
        }
        purrr::flatten(current_record)
      }
    ) %>%
    dplyr::select(
      "taxon", "guid", "mbNumber", "taxonomicLevel", "trophicMode",
      "guild", "confidenceRanking", "growthForm", "trait", "notes",
      "citationSource"
    )
}

#' Assign Guilds to Organisms Based on Taxonomic Classification
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' The original function and documentation was written by Brendan Furneaux
#' in the [FUNGuildR](https://github.com/brendanf/FUNGuildR/) package.
#'
#'
#' These functions have identical behavior if supplied with a database; however
#' they download the database corresponding to their name by default.
#'
#' Taxa present in the database are matched to the taxa present in the supplied
#' `otu_table` by exact name.
#' In the case of multiple matches, the lowest (most specific) rank is chosen.
#' No attempt is made to check or correct the classification in
#' `otu_table$Taxonomy`.
#'
#' @param otu_table A `data.frame` with a `character`
#' column named "`Taxonomy`" (or another name as specified in
#' `tax_col`), as well as any other columns.
#' Each entry in "`otu_table$Taxonomy`" should be a comma-, colon-,
#' underscore-, or semicolon-delimited classification of an organism.
#' Rank indicators as given by Sintax ("`k:`", "`p:`"...) or Unite ("`k__`,
#' "`p__`", ...) are also allowed.
#' A `character` vector, representing only the taxonomic classification,
#' is also accepted.
#' @param tax_col A `character` string, optionally giving an alternate
#' column name in `otu_table` to use instead of `otu_table$Taxonomy`.
#'
#' @param db_funguild A `data.frame` representing the FUNGuild as returned by
#' [get_funguild_db()]
#' If not supplied, the default database will be downloaded.
#'
#' @return A [`tibble::tibble`] containing all columns of
#' `otu_table`, plus relevant columns of information from the FUNGuild
#' @export
#'
#' @references Nguyen NH, Song Z, Bates ST, Branco S, Tedersoo L, Menke J,
#' Schilling JS, Kennedy PG. 2016. *FUNGuild: An open annotation tool for
#' parsing fungal community datasets by ecological guild*. Fungal Ecology
#' 20:241-248.
#' @author Brendan Furneaux (orcid: [0000-0003-3522-7363](https://orcid.org/0000-0003-3522-7363)),
#' modified by Adrien Taudière
funguild_assign <- function(
        otu_table, db_funguild = get_funguild_db(),
        tax_col = "Taxonomy") {
  if (is.character(otu_table)) {
    otu_table <- tibble::tibble(otu_table)
    names(otu_table) <- tax_col
  }
  if (!is.data.frame(otu_table)) {
    stop(paste0("otu_table must be a dataframe not a ", class(otu_table)))
  }
  if (!tax_col %in% colnames(otu_table)) {
    stop(paste0("The tax_col args ", tax_col, " is not present in the colnames
    of the otu_table dataframe."))
  }

  make_taxkey <- function(x) {
    out <- gsub("\\b[kpcofgs](:|__)", "", x)
    out <- gsub("[_ ;,:]", "@", out)
    paste0("@", out, "@")
  }

  otu_table$taxkey <- make_taxkey(otu_table[[tax_col]])
  all_taxkey <- unique(otu_table$taxkey) %>% na.omit()
  `.` <- taxon <- taxkey <- searchkey <- taxonomicLevel <- NULL # to pass R CMD check
  db_funguild <- dplyr::mutate(
    db_funguild,
    searchkey = paste0("@", stringr::str_replace(taxon, "[ _]", "@"), "@")
  )
  dplyr::select(db_funguild, taxonomicLevel, searchkey) %>%
    dplyr::mutate(
      taxkey = purrr::map(searchkey, stringr::str_subset, string = all_taxkey)
    ) %>%
    tidyr::unnest(cols = taxkey) %>%
    dplyr::group_by(taxkey) %>%
    dplyr::arrange(dplyr::desc(taxonomicLevel)) %>%
    dplyr::summarize_at("searchkey", dplyr::first) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_all(as.character) %>%
    dplyr::left_join(otu_table, ., by = "taxkey") %>%
    dplyr::left_join(db_funguild, by = "searchkey", suffix = c("", ".funguild")) %>%
    dplyr::select(-taxkey, -searchkey)
}
